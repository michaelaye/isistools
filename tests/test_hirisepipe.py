"""Tests for the hirisepipe HiRISE calibration pipeline.

Compares Python hical output against ISIS reference cubes.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import xarray as xr

# Test data paths
_REF_DIR = Path(__file__).parent / "data" / "hirise_reference"
_REF_RAW = _REF_DIR / "raw.cub"
_REF_CAL = _REF_DIR / "cal.cub"
_HIRISE_EDR = Path.home() / "Dropbox" / "data" / "hirise" / "ESP_021491_0950_RED4_0.IMG"

pytestmark = pytest.mark.skipif(
    not _REF_RAW.exists(),
    reason="HiRISE test data not available",
)


def _read_isis_cube(path: Path) -> np.ndarray:
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Dataset has no geotransform")
        da = xr.open_dataarray(str(path), engine="rasterio")
    return da.values.squeeze()


class TestHiCalParams:
    def test_from_cube(self):
        from isistools.hirisepipe.hical import HiCalParams

        params = HiCalParams.from_cube(_REF_RAW)
        assert params.ccd == 4
        assert params.ccd_name == "RED4"
        assert params.channel == 0
        assert params.tdi == 128
        assert params.binning == 1
        assert params.n_samples == 1024
        assert params.n_lines == 35000
        assert params.filter_name == "RED"


class TestHiCal:
    @pytest.fixture(scope="class")
    def calibrated(self):
        from isistools.hirisepipe.hical import hical

        return hical(_REF_RAW, units="DN")

    def test_shape(self, calibrated):
        assert calibrated.shape == (35000, 1024)

    def test_dtype(self, calibrated):
        assert calibrated.dtype == np.float32

    @pytest.mark.skipif(not _REF_CAL.exists(), reason="ISIS hical reference not available")
    def test_matches_isis_mean(self, calibrated):
        """Mean must match ISIS within 0.1%."""
        isis_cal = _read_isis_cube(_REF_CAL)
        py_mean = np.nanmean(calibrated)
        isis_mean = np.nanmean(isis_cal)
        rel_err = abs(py_mean - isis_mean) / isis_mean
        assert rel_err < 0.001, f"Mean rel error {rel_err:.6f} > 0.1%"

    @pytest.mark.skipif(not _REF_CAL.exists(), reason="ISIS hical reference not available")
    def test_matches_isis_std(self, calibrated):
        """Std dev must match ISIS within 0.5%."""
        isis_cal = _read_isis_cube(_REF_CAL)
        py_std = np.nanstd(calibrated)
        isis_std = np.nanstd(isis_cal)
        rel_err = abs(py_std - isis_std) / isis_std
        assert rel_err < 0.005, f"StdDev rel error {rel_err:.6f} > 0.5%"

    @pytest.mark.skipif(not _REF_CAL.exists(), reason="ISIS hical reference not available")
    def test_max_relative_error(self, calibrated):
        """Max per-pixel relative error must be < 0.5%."""
        isis_cal = _read_isis_cube(_REF_CAL)
        valid = np.isfinite(calibrated) & np.isfinite(isis_cal)
        diff = np.abs(calibrated[valid].astype(np.float64) - isis_cal[valid].astype(np.float64))
        rel = diff / np.abs(isis_cal[valid].astype(np.float64))
        assert rel.max() < 0.005, f"Max rel error {rel.max():.6f} > 0.5%"

    def test_stats_sanity(self, calibrated):
        """Calibrated RED4 should have reasonable DN values."""
        valid = calibrated[np.isfinite(calibrated)]
        assert 1000 < valid.min() < 1200
        assert 2500 < valid.max() < 3200
        assert 1400 < valid.mean() < 1800


class TestStitch:
    def test_stitch_basic(self):
        from isistools.hirisepipe.stitch import stitch_channels

        ch0 = np.full((100, 50), 100.0, dtype=np.float32)
        ch1 = np.full((100, 50), 200.0, dtype=np.float32)
        result = stitch_channels(ch0, ch1, balance=False)
        assert result.shape == (100, 100)
        # ch1 on left, ch0 on right
        assert result[0, 0] == 200.0  # ch1
        assert result[0, 50] == 100.0  # ch0

    def test_stitch_balance(self):
        from isistools.hirisepipe.stitch import stitch_channels

        ch0 = np.full((100, 50), 100.0, dtype=np.float32)
        ch1 = np.full((100, 50), 200.0, dtype=np.float32)
        result = stitch_channels(ch0, ch1, balance=True, truth_channel=0)
        # ch1 should be scaled to match ch0's level
        # Since seam averages: ch0_seam=100, ch1_seam=200, coeff=0.5
        assert result.shape == (100, 100)
        np.testing.assert_allclose(result[0, 0], 100.0, atol=0.01)  # ch1 scaled
        assert result[0, 50] == 100.0  # ch0 unchanged


class TestCubenorm:
    def test_cubenorm_preserves_mean(self):
        from isistools.hirisepipe.cubenorm import cubenorm

        rng = np.random.default_rng(42)
        image = rng.normal(1000, 50, (200, 50)).astype(np.float32)
        # Add column variation
        image *= np.linspace(0.8, 1.2, 50)[np.newaxis, :]

        result = cubenorm(image, preserve=True)
        np.testing.assert_allclose(
            np.nanmean(result), np.nanmean(image), rtol=0.01
        )

    def test_cubenorm_reduces_column_variation(self):
        from isistools.hirisepipe.cubenorm import cubenorm

        rng = np.random.default_rng(42)
        image = rng.normal(1000, 10, (200, 50)).astype(np.float32)
        # Add strong column variation
        col_factors = np.linspace(0.5, 1.5, 50)
        image *= col_factors[np.newaxis, :]

        result = cubenorm(image)

        # Column std should be much smaller after normalization
        col_means_before = np.nanmean(image, axis=0)
        col_means_after = np.nanmean(result, axis=0)
        assert col_means_after.std() < col_means_before.std() * 0.1


class TestIngestEDR:
    """Test direct PDS EDR ingestion (no ISIS needed)."""

    @pytest.mark.skipif(
        not _HIRISE_EDR.exists(),
        reason="HiRISE EDR not available",
    )
    def test_ingest_matches_isis(self):
        """EDR reader must produce pixel-exact match with hi2isis output."""
        from isistools.hirisepipe.ingest import ingest_hirise_edr

        edr = ingest_hirise_edr(_HIRISE_EDR)
        isis_raw = _read_isis_cube(_REF_RAW)
        np.testing.assert_array_equal(edr.image.astype(np.float32), isis_raw)

    @pytest.mark.skipif(
        not _HIRISE_EDR.exists(),
        reason="HiRISE EDR not available",
    )
    def test_ingest_metadata(self):
        from isistools.hirisepipe.ingest import ingest_hirise_edr

        edr = ingest_hirise_edr(_HIRISE_EDR)
        assert edr.cpmm_number == 5
        assert edr.channel_number == 0
        assert edr.tdi == 128
        assert edr.binning == 1
        assert edr.n_samples == 1024
        assert edr.n_lines == 35000
        assert edr.filter_name == "RED"

    @pytest.mark.skipif(
        not _HIRISE_EDR.exists() or not _REF_CAL.exists(),
        reason="HiRISE test data not available",
    )
    def test_hical_from_edr_matches_isis(self):
        """hical_from_edr must match ISIS hical within 0.5%."""
        from isistools.hirisepipe.hical import hical_from_edr

        cal = hical_from_edr(_HIRISE_EDR, units="DN")
        isis_cal = _read_isis_cube(_REF_CAL)

        valid = np.isfinite(cal) & np.isfinite(isis_cal)
        rel = np.abs(cal[valid].astype(np.float64) - isis_cal[valid].astype(np.float64))
        rel /= np.abs(isis_cal[valid].astype(np.float64))
        assert rel.max() < 0.005, f"Max rel error {rel.max():.6f} > 0.5%"
