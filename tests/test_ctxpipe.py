"""Tests for the ctxpipe CTX calibration pipeline.

Compares Python pipeline output against ISIS reference cubes at each stage.
The reference cubes are produced by running the ISIS CTX pipeline:
    mroctx2isis -> ctxcal (IOF=false) -> ctxevenodd
on the test EDR B04_011267_0983_XN_81S063W.IMG.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from isistools.ctxpipe.calibrate import calibrate
from isistools.ctxpipe.evenodd import correct_evenodd
from isistools.ctxpipe.ingest import CTXMetadata, ingest_ctx_edr

# Test data paths
_DATA_DIR = Path(__file__).parent / "data"
_CTX_EDR = Path.home() / "Dropbox" / "data" / "ctx" / "B04_011267_0983_XN_81S063W.IMG"
_REF_DIR = _DATA_DIR / "ctx_reference"
_REF_RAW = _REF_DIR / "raw.cub"
_REF_CAL_IOF = _REF_DIR / "cal_iof.cub"
_REF_CAL = _REF_DIR / "cal.cub"
_REF_EVENODD = _REF_DIR / "eveNodd.cub"

# Skip all tests if EDR or reference cubes are missing
pytestmark = pytest.mark.skipif(
    not _CTX_EDR.exists() or not _REF_RAW.exists(),
    reason="CTX test data not available",
)


def _read_isis_cube(path: Path) -> np.ndarray:
    """Read an ISIS cube via rasterio and return the 2D array."""
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Dataset has no geotransform")
        da = xr.open_dataarray(str(path), engine="rasterio")
    return da.values.squeeze()


@pytest.fixture(scope="module")
def ingested():
    """Ingest the CTX EDR once for all tests."""
    return ingest_ctx_edr(_CTX_EDR)


@pytest.fixture(scope="module")
def calibrated(ingested):
    """Calibrate the ingested image once for all tests."""
    image, meta = ingested
    return calibrate(image, meta, iof=False)


@pytest.fixture(scope="module")
def evenodd_corrected(calibrated, ingested):
    """Apply even/odd correction once for all tests."""
    _, meta = ingested
    return correct_evenodd(calibrated, spatial_summing=meta.spatial_summing)


class TestIngest:
    """Test CTX EDR ingestion against ISIS mroctx2isis."""

    def test_shape(self, ingested):
        image, meta = ingested
        assert image.shape == (11264, 5000)

    def test_pixel_exact_match(self, ingested):
        """Python ingestion must match ISIS mroctx2isis pixel-for-pixel."""
        image, _ = ingested
        isis_raw = _read_isis_cube(_REF_RAW)
        # Compare as int16 since ISIS stores as SignedWord
        np.testing.assert_array_equal(
            image.astype(np.int16),
            isis_raw.astype(np.int16),
        )

    def test_metadata(self, ingested):
        _, meta = ingested
        assert meta.product_id == "B04_011267_0983_XN_81S063W"
        assert meta.spatial_summing == 1
        assert meta.sample_first_pixel == 0
        assert meta.line_exposure_duration == 1.877
        assert meta.target_name.upper() == "MARS"

    def test_dark_pixels_shape(self, ingested):
        _, meta = ingested
        assert meta.dark_pixels.shape == (11264, 24)

    def test_dark_pixels_match_isis(self, ingested):
        """Dark pixels must match the ISIS 'Ctx Prefix Dark Pixels' table."""
        _, meta = ingested
        # Read the dark pixel table from the ISIS cube
        with open(_REF_RAW, "rb") as f:
            f.seek(112706113 - 1)  # ISIS 1-based StartByte
            isis_dark = np.frombuffer(
                f.read(1081344), dtype="<i4"
            ).reshape(11264, 24)
        np.testing.assert_array_equal(meta.dark_pixels.astype(np.int32), isis_dark)

    def test_no_nans(self, ingested):
        """This particular EDR has no DN=0 pixels, so no NaNs expected."""
        image, _ = ingested
        assert not np.isnan(image).any()


class TestCalibrate:
    """Test CTX calibration against ISIS ctxcal (IOF=false)."""

    def test_shape(self, calibrated):
        assert calibrated.shape == (11264, 5000)

    @pytest.mark.skipif(
        not _REF_CAL.exists(),
        reason="ISIS ctxcal reference cube not available",
    )
    def test_matches_isis_ctxcal(self, calibrated):
        """Calibrated output must match ISIS ctxcal to float32 precision."""
        isis_cal = _read_isis_cube(_REF_CAL)
        # Allow float32 machine epsilon tolerance
        np.testing.assert_allclose(
            calibrated, isis_cal, rtol=1e-6, atol=1e-5
        )

    @pytest.mark.skipif(
        not _REF_CAL.exists(),
        reason="ISIS ctxcal reference cube not available",
    )
    def test_max_relative_error(self, calibrated):
        """Max relative error must be < 1e-6 (float32 precision)."""
        isis_cal = _read_isis_cube(_REF_CAL)
        rel_err = np.abs(calibrated - isis_cal) / np.abs(isis_cal)
        assert np.nanmax(rel_err) < 1e-6

    def test_stats(self, calibrated):
        """Sanity check on calibrated statistics."""
        valid = calibrated[np.isfinite(calibrated)]
        assert 14.0 < valid.min() < 15.0
        assert 56.0 < valid.max() < 58.0
        assert 27.0 < valid.mean() < 29.0


class TestEvenOdd:
    """Test even/odd correction against ISIS ctxevenodd."""

    def test_shape(self, evenodd_corrected):
        assert evenodd_corrected.shape == (11264, 5000)

    @pytest.mark.skipif(
        not _REF_EVENODD.exists(),
        reason="ISIS ctxevenodd reference cube not available",
    )
    def test_matches_isis_ctxevenodd(self, evenodd_corrected):
        """Even/odd output must match ISIS ctxevenodd to float32 precision."""
        isis_eo = _read_isis_cube(_REF_EVENODD)
        np.testing.assert_allclose(
            evenodd_corrected, isis_eo, rtol=1e-6, atol=1e-5
        )

    def test_reduces_striping(self, calibrated, evenodd_corrected):
        """Even/odd correction should reduce the difference between
        odd and even column means."""
        # Before correction
        cal_odd_mean = np.nanmean(calibrated[:, 0::2])
        cal_even_mean = np.nanmean(calibrated[:, 1::2])
        before_diff = abs(cal_odd_mean - cal_even_mean)

        # After correction
        eo_odd_mean = np.nanmean(evenodd_corrected[:, 0::2])
        eo_even_mean = np.nanmean(evenodd_corrected[:, 1::2])
        after_diff = abs(eo_odd_mean - eo_even_mean)

        assert after_diff < before_diff

    def test_preserves_mean(self, calibrated, evenodd_corrected):
        """Overall mean should be preserved by even/odd correction."""
        np.testing.assert_allclose(
            np.nanmean(calibrated),
            np.nanmean(evenodd_corrected),
            rtol=1e-5,
        )

    def test_skipped_for_summing_2(self):
        """Even/odd correction should be a no-op for summing > 1."""
        dummy = np.random.rand(10, 20).astype(np.float32)
        result = correct_evenodd(dummy, spatial_summing=2)
        np.testing.assert_array_equal(result, dummy)


class TestPipeline:
    """Test the full pipeline end-to-end."""

    def test_full_pipeline(self):
        from isistools.ctxpipe.pipeline import ctx_calibrate

        image, meta = ctx_calibrate(_CTX_EDR, iof=False)
        assert image.shape == (11264, 5000)
        assert meta.product_id == "B04_011267_0983_XN_81S063W"

        valid = image[np.isfinite(image)]
        assert 14.0 < valid.min() < 16.0
        assert 55.0 < valid.max() < 58.0

    def test_iof_requires_distance(self):
        from isistools.ctxpipe.calibrate import calibrate

        dummy_img = np.ones((10, 10), dtype=np.float32)
        dummy_meta = CTXMetadata(
            product_id="test",
            spacecraft_clock_count="",
            start_time="",
            target_name="Mars",
            mission_phase_name="",
            spatial_summing=1,
            sample_first_pixel=0,
            line_exposure_duration=1.0,
            focal_plane_temperature=290.0,
            offset_mode_id="",
            sample_bit_mode_id="SQROOT",
            lines=10,
            line_samples=10,
            data_set_id="",
            orbit_number=0,
        )
        dummy_meta.dark_pixels = np.zeros((10, 24), dtype=np.int16)

        with pytest.raises(ValueError, match="sun_distance_km is required"):
            calibrate(dummy_img, dummy_meta, iof=True)


class TestIoF:
    """Test I/F calibration against ISIS ctxcal IOF=true."""

    @pytest.mark.skipif(
        not _REF_CAL_IOF.exists(),
        reason="ISIS ctxcal IOF reference cube not available",
    )
    def test_iof_matches_isis(self, ingested):
        """I/F calibration must match ISIS to float32 precision."""
        image, meta = ingested
        cal_iof = calibrate(image, meta, iof=True, sun_distance_km=220104537.28)
        isis_iof = _read_isis_cube(_REF_CAL_IOF)
        np.testing.assert_allclose(cal_iof, isis_iof, rtol=1e-6, atol=1e-8)

    @pytest.mark.skipif(
        not _REF_RAW.exists(),
        reason="CTX test data not available",
    )
    def test_sun_distance_from_cube(self):
        """Auto Sun distance computation via spiceypy."""
        pytest.importorskip("spiceypy")
        from isistools.ctxpipe.spice_utils import sun_distance_from_cube

        dist = sun_distance_from_cube(_REF_RAW)
        assert 2.0e8 < dist < 2.5e8  # ~1.3-1.7 AU in km
