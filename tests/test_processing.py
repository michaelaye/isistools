"""Tests for the csm2map processing pipeline.

Unit tests for grid, resample, and special pixel handling run without
CSM/ale dependencies. Integration tests that need real data are marked
with @pytest.mark.slow.
"""

import numpy as np
import pytest


def test_grid_from_params():
    """Test basic grid construction from explicit parameters."""
    from isistools.processing.grid import grid_from_params

    grid = grid_from_params(
        crs="+proj=eqc +lat_ts=0 +lon_0=0 +a=3396190 +b=3376200 +units=m +no_defs +type=crs",
        resolution=6.0,
        lat_min=-5.0,
        lat_max=-3.0,
        lon_min=130.0,
        lon_max=132.0,
    )
    assert grid.width > 0
    assert grid.height > 0
    assert grid.resolution == 6.0
    assert grid.lat_min == -5.0
    assert grid.lon_max == 132.0


def test_resample_identity():
    """Test resampling with an identity coordinate map."""
    from isistools.processing.resample import Interpolation, resample
    from isistools.processing.transform import CoordinateMap

    # Create a simple 10x10 test image
    input_data = np.arange(100, dtype=np.float32).reshape(10, 10)

    # Identity mapping: output pixel (r, c) maps to input pixel (r, c)
    rows, cols = np.meshgrid(np.arange(10), np.arange(10), indexing="ij")
    coord_map = CoordinateMap(
        input_lines=rows.astype(np.float64),
        input_samples=cols.astype(np.float64),
        valid=np.ones((10, 10), dtype=bool),
    )

    result = resample(input_data, coord_map, interpolation=Interpolation.NEAREST)
    np.testing.assert_array_almost_equal(result, input_data, decimal=5)


def test_resample_shift():
    """Test resampling with a half-pixel shift."""
    from isistools.processing.resample import Interpolation, resample
    from isistools.processing.transform import CoordinateMap

    input_data = np.ones((10, 10), dtype=np.float32) * 42.0

    rows, cols = np.meshgrid(np.arange(10), np.arange(10), indexing="ij")
    coord_map = CoordinateMap(
        input_lines=rows.astype(np.float64) + 0.5,
        input_samples=cols.astype(np.float64) + 0.5,
        valid=np.ones((10, 10), dtype=bool),
    )

    result = resample(input_data, coord_map, interpolation=Interpolation.BILINEAR)
    # Constant image shifted should still be ~42 everywhere (except edges)
    interior = result[1:-1, 1:-1]
    np.testing.assert_array_almost_equal(interior, 42.0, decimal=3)


def test_isis_special_pixel_masking():
    """Test that ISIS special pixels are converted to NaN."""
    from isistools.io.cubes import ISIS_NULL, _mask_special_pixels

    data = np.array([1.0, 2.0, ISIS_NULL, 3.0], dtype=np.float32)
    _mask_special_pixels(data)
    assert np.isnan(data[2])
    assert data[0] == 1.0


def test_mapping_to_crs():
    """Test that mapping_to_crs returns a valid pyproj CRS."""
    from isistools.geo.projections import mapping_to_crs

    mapping = {
        "ProjectionName": "Equirectangular",
        "EquatorialRadius": 3396190.0,
        "PolarRadius": 3376200.0,
        "CenterLongitude": 0.0,
        "CenterLatitude": 0.0,
    }
    crs = mapping_to_crs(mapping)
    # Should be a valid CRS with equirectangular projection
    proj4 = crs.to_proj4()
    assert "+proj=eqc" in proj4


@pytest.mark.slow
def test_ctx_end_to_end():
    """End-to-end test with real CTX data. Requires ISIS + CSM deps + data."""
    pytest.skip("Requires CTX test data -- run manually")
