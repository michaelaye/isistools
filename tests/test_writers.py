"""Nodata handling in the csm2map GeoTIFF writer.

Regression cover for the case where 0.0 was declared as nodata. 0.0 is a legal
pixel value, so every genuinely-zero pixel was indistinguishable from empty: on
a 20-image MC20W set, 0.25 % of the pixels cam2map marked valid were exactly
0.0 and would have been dropped by any reader honouring the declared nodata.
"""

from __future__ import annotations

import numpy as np
import pytest
import rasterio
from affine import Affine
from rasterio.crs import CRS as RioCRS

from isistools.csm2map.grid import OutputGrid
from isistools.csm2map.writers import write_geotiff


@pytest.fixture
def grid():
    """A 4x4 grid; the projection is irrelevant to nodata handling."""
    return OutputGrid(
        crs=RioCRS.from_epsg(4326),
        transform=Affine(5.0, 0.0, 0.0, 0.0, -5.0, 0.0),
        width=4,
        height=4,
        x_min=0.0,
        x_max=20.0,
        y_min=-20.0,
        y_max=0.0,
        resolution=5.0,
        lat_min=-1.0,
        lat_max=0.0,
        lon_min=0.0,
        lon_max=1.0,
    )


@pytest.fixture
def data():
    """Two pixels that matter: a real 0.0, and a NaN standing for no data."""
    a = np.arange(16, dtype=np.float32).reshape(4, 4) / 100.0
    a[0, 0] = 0.0  # legal value, must survive as data
    a[3, 3] = np.nan  # no data
    return a


def test_default_nodata_is_nan(tmp_path, grid, data):
    out = write_geotiff(tmp_path / "o.tif", data, grid)
    with rasterio.open(out) as d:
        assert np.isnan(d.nodata)


def test_a_real_zero_is_not_treated_as_nodata(tmp_path, grid, data):
    """The actual regression: 0.0 must read back as data, not as empty."""
    out = write_geotiff(tmp_path / "o.tif", data, grid)
    with rasterio.open(out) as d:
        masked = d.read(1, masked=True)
    assert masked[0, 0] == 0.0
    assert not masked.mask[0, 0], "a legal 0.0 pixel was masked as nodata"


def test_nan_pixels_are_the_masked_ones(tmp_path, grid, data):
    out = write_geotiff(tmp_path / "o.tif", data, grid)
    with rasterio.open(out) as d:
        raw = d.read(1)
        masked = d.read(1, masked=True)
    assert np.array_equal(masked.mask, np.isnan(raw))
    assert masked.mask[3, 3]
    assert masked.mask.sum() == 1


def test_an_explicit_nodata_is_still_honoured(tmp_path, grid, data):
    """--nodata exists for consumers that cannot handle NaN, so the
    substitution path must still work."""
    out = write_geotiff(tmp_path / "o.tif", data, grid, nodata=-9999.0)
    with rasterio.open(out) as d:
        raw = d.read(1)
        assert d.nodata == -9999.0
    assert raw[3, 3] == -9999.0
    assert not np.isnan(raw).any(), "NaN should have been substituted"
    assert raw[0, 0] == 0.0, "a real zero must survive a non-NaN nodata too"
