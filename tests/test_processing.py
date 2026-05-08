"""Tests for the csm2map processing pipeline.

Unit tests for grid, resample, and special pixel handling run without
CSM/ale dependencies. Integration tests that need real data are marked
with @pytest.mark.slow.
"""

import numpy as np
import pytest


def test_grid_from_params():
    """Test basic grid construction from explicit parameters."""
    from isistools.csm2map.grid import grid_from_params

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
    from isistools.csm2map.resample import Interpolation, resample
    from isistools.csm2map.transform import CoordinateMap

    # Create a simple 10x10 test image
    input_data = np.arange(100, dtype=np.float32).reshape(10, 10)

    # Identity mapping: output pixel (r, c) maps to input pixel (r, c)
    rows, cols = np.meshgrid(np.arange(10), np.arange(10), indexing="ij")
    coords = np.stack(
        (rows.astype(np.float64), cols.astype(np.float64)),
        axis=0,
    )
    coord_map = CoordinateMap(coords=coords, valid=np.ones((10, 10), dtype=bool))

    result = resample(input_data, coord_map, interpolation=Interpolation.NEAREST)
    np.testing.assert_array_almost_equal(result, input_data, decimal=5)


def test_resample_shift():
    """Test resampling with a half-pixel shift."""
    from isistools.csm2map.resample import Interpolation, resample
    from isistools.csm2map.transform import CoordinateMap

    input_data = np.ones((10, 10), dtype=np.float32) * 42.0

    rows, cols = np.meshgrid(np.arange(10), np.arange(10), indexing="ij")
    coords = np.stack(
        (rows.astype(np.float64) + 0.5, cols.astype(np.float64) + 0.5),
        axis=0,
    )
    coord_map = CoordinateMap(coords=coords, valid=np.ones((10, 10), dtype=bool))

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
    from isistools.csm2map.projections import mapping_to_crs

    mapping = {
        "ProjectionName": "Equirectangular",
        "EquatorialRadius": 3396190.0,
        "PolarRadius": 3376200.0,
        "CenterLongitude": 0.0,
        "CenterLatitude": 0.0,
    }
    crs = mapping_to_crs(mapping)
    # Should be a valid CRS with equirectangular projection
    assert crs.coordinate_operation is not None
    assert "Cylindrical" in crs.coordinate_operation.method_name


def test_coordinatemap_view_semantics():
    """input_lines / input_samples must be views over coords, not copies.

    Catches future refactors that accidentally drop the view semantics
    of the consolidated `(2, h, w)` storage. If a copy creeps in, the
    per-stripe slice in `_resample_band` reintroduces the
    multi-GB transient that L1 was specifically designed to avoid.
    """
    from isistools.csm2map.transform import CoordinateMap

    coords = np.zeros((2, 4, 5), dtype=np.float32)
    coord_map = CoordinateMap(coords=coords, valid=np.ones((4, 5), dtype=bool))

    assert np.shares_memory(coord_map.input_lines, coord_map.coords)
    assert np.shares_memory(coord_map.input_samples, coord_map.coords)

    # Mutating through the property must show up in the underlying buffer.
    coord_map.input_lines[0, 0] = 7.0
    coord_map.input_samples[1, 1] = 9.0
    assert coords[0, 0, 0] == 7.0
    assert coords[1, 1, 1] == 9.0
    assert coord_map.shape == (4, 5)


def test_compute_transform_dense_validity_mask(monkeypatch):
    """Dense transform validity mask: NaN, in-bounds, out-of-bounds.

    Mocks `ground_to_image_batch` to return a deterministic 3x3 grid
    so we can assert the chained `&=` validity logic produces the
    correct mask without needing a real CSM model.
    """
    from isistools.csm2map import transform as transform_mod
    from isistools.csm2map.grid import grid_from_params

    grid = grid_from_params(
        crs="+proj=eqc +lat_ts=0 +lon_0=0 +a=3396190 +b=3376200 +units=m +no_defs +type=crs",
        resolution=10000.0,
        lat_min=-1.0,
        lat_max=1.0,
        lon_min=-1.0,
        lon_max=1.0,
    )

    def fake_ground_to_image(model, lat, lon, radii):
        h, w = lat.shape
        lines = np.full((h, w), 100.0)
        samps = np.full((h, w), 100.0)
        # NaN at (0, 0) — fails isfinite
        lines[0, 0] = np.nan
        # Out-of-bounds line at (h-1, w-1)
        lines[-1, -1] = 999_999.0
        # Out-of-bounds sample at (0, w-1)
        samps[0, -1] = -10.0
        return lines, samps

    monkeypatch.setattr(transform_mod, "ground_to_image_batch", fake_ground_to_image)

    coord_map = transform_mod.compute_transform_dense(
        model=None,
        grid=grid,
        surface_radius=3396190.0,
        input_n_lines=500,
        input_n_samples=500,
    )
    assert not coord_map.valid[0, 0]  # NaN
    assert not coord_map.valid[-1, -1]  # OOB line
    assert not coord_map.valid[0, -1]  # OOB sample
    # Interior pixels should be valid (line=samp=100, in bounds [0, 500))
    assert coord_map.valid[1, 0]


def test_compute_transform_coarse_validity_mask(monkeypatch):
    """Coarse transform validity mask uses chained `&=` correctly.

    Uses a small grid where step >= grid size, so the coarse path
    degenerates to a 2x2 coarse grid + bilinear upsample to the full
    output, exercising the same bool chain as the dense path.
    """
    from isistools.csm2map import transform as transform_mod
    from isistools.csm2map.grid import grid_from_params

    grid = grid_from_params(
        crs="+proj=eqc +lat_ts=0 +lon_0=0 +a=3396190 +b=3376200 +units=m +no_defs +type=crs",
        resolution=10000.0,
        lat_min=-1.0,
        lat_max=1.0,
        lon_min=-1.0,
        lon_max=1.0,
    )

    def fake_ground_to_image(model, lat, lon, radii):
        # Returns the same shape as input lat/lon — coarse grid
        return np.full(lat.shape, 100.0), np.full(lat.shape, 100.0)

    monkeypatch.setattr(transform_mod, "ground_to_image_batch", fake_ground_to_image)

    coord_map = transform_mod.compute_transform_coarse(
        model=None,
        grid=grid,
        surface_radius=3396190.0,
        step=64,  # grid is small; coarse grid degenerates to corners
        input_n_lines=500,
        input_n_samples=500,
    )
    # Constant 100,100 input is in-bounds [0, 500), so all valid
    assert coord_map.valid.all()
    assert coord_map.shape == (grid.height, grid.width)
    # And the consolidated coords storage was used
    assert coord_map.coords.shape == (2, grid.height, grid.width)


def test_bilinear_upsample_window_matches_full():
    """Windowed bilinear upsample is bit-identical to slicing a full upsample.

    This is the key correctness invariant for the tiled path: tile
    boundaries must be pixel-exact, so adjacent tiles agree without any
    overlap-and-trim. Equivalent to asserting the L2 plan's "build
    coarse grid once globally, windowed upsample per tile" yields the
    same pixels as a single global upsample.
    """
    from isistools.csm2map.transform import (
        _bilinear_upsample_pair,
        _bilinear_upsample_pair_window,
    )

    rng = np.random.default_rng(42)
    nrc, ncc = 5, 7
    coarse_a = rng.standard_normal((nrc, ncc)).astype(np.float32)
    coarse_b = rng.standard_normal((nrc, ncc)).astype(np.float32)
    full_h, full_w = 73, 91

    full = _bilinear_upsample_pair(coarse_a, coarse_b, full_h, full_w)

    # Tile the (full_h, full_w) output into 4 windows of varying sizes,
    # including the bottom-right edge tile that's not a multiple of
    # tile_size, and assert each matches the corresponding slice of full.
    tiles = [
        (0, 0, 40, 50),
        (0, 50, 40, 41),
        (40, 0, 33, 50),
        (40, 50, 33, 41),
    ]
    for row0, col0, h_t, w_t in tiles:
        windowed = _bilinear_upsample_pair_window(
            coarse_a, coarse_b, full_h, full_w, row0, col0, h_t, w_t
        )
        expected = full[:, row0 : row0 + h_t, col0 : col0 + w_t]
        np.testing.assert_array_equal(windowed, expected)


def test_coordinate_map_for_window_matches_full(monkeypatch):
    """coordinate_map_for_window slices match a full compute_transform_coarse.

    Builds a small synthetic grid, runs both paths through a
    monkeypatched ground_to_image_batch, and asserts that for any
    window the tiled path produces the same coords + valid mask as the
    corresponding slice of the global CoordinateMap.
    """
    from isistools.csm2map import transform as transform_mod
    from isistools.csm2map.grid import grid_from_params

    grid = grid_from_params(
        crs="+proj=eqc +lat_ts=0 +lon_0=0 +a=3396190 +b=3376200 +units=m +no_defs +type=crs",
        resolution=10000.0,
        lat_min=-2.0,
        lat_max=2.0,
        lon_min=-2.0,
        lon_max=2.0,
    )

    def fake_ground_to_image(model, lat, lon, radii):
        # Smoothly varying lines/samples driven by lat/lon so the
        # bilinear upsample exercises non-trivial gradients. Sprinkle
        # one NaN and one OOB value so the validity mask is non-trivial.
        lines = (300.0 + 50.0 * lat / 0.04).astype(np.float32)
        samps = (300.0 + 50.0 * lon / 0.04).astype(np.float32)
        # Force one NaN at index (0, 0)
        if lines.size > 0:
            lines.flat[0] = np.nan
        # Force one OOB sample (large) somewhere predictable
        if samps.size > 1:
            samps.flat[-1] = 1_000_000.0
        return lines, samps

    monkeypatch.setattr(transform_mod, "ground_to_image_batch", fake_ground_to_image)

    full = transform_mod.compute_transform_coarse(
        model=None,
        grid=grid,
        surface_radius=3396190.0,
        step=8,
        input_n_lines=500,
        input_n_samples=500,
    )
    state = transform_mod.compute_coarse_state(
        model=None,
        grid=grid,
        surface_radius=3396190.0,
        step=8,
        input_n_lines=500,
        input_n_samples=500,
    )

    h, w = grid.height, grid.width
    # A window that covers the full grid in two halves, plus an interior
    # block, exercises tile boundaries and edge handling.
    h_half = h // 2
    w_half = w // 2
    windows = [
        (0, 0, h_half, w_half),
        (0, w_half, h_half, w - w_half),
        (h_half, 0, h - h_half, w_half),
        (h_half, w_half, h - h_half, w - w_half),
    ]
    for row0, col0, h_t, w_t in windows:
        win_map = transform_mod.coordinate_map_for_window(state, row0, col0, h_t, w_t)
        np.testing.assert_array_equal(
            win_map.coords, full.coords[:, row0 : row0 + h_t, col0 : col0 + w_t]
        )
        np.testing.assert_array_equal(
            win_map.valid, full.valid[row0 : row0 + h_t, col0 : col0 + w_t]
        )


def test_resolve_tile_size_none_forces_batch():
    """``"none"`` / 0 / None all return ``None`` → batch path."""
    from isistools.csm2map.tiled import resolve_tile_size

    for spec in (None, 0, "none", "None"):
        assert (
            resolve_tile_size(spec, grid_h=10000, grid_w=10000, persistent_rss_bytes=10**9)
            is None
        )


def test_resolve_tile_size_int_passthrough():
    """A positive int spec is returned verbatim, no auto fall-through."""
    from isistools.csm2map.tiled import resolve_tile_size

    assert (
        resolve_tile_size(4096, grid_h=100000, grid_w=100000, persistent_rss_bytes=10**9)
        == 4096
    )
    # Even a tile larger than the output is honored — caller asked explicitly.
    assert (
        resolve_tile_size(20000, grid_h=10, grid_w=10, persistent_rss_bytes=10**9) == 20000
    )


def test_resolve_tile_size_auto_falls_through_to_batch_when_tile_covers_output():
    """``"auto"`` returns ``None`` when the chosen tile covers the full output."""
    from isistools.csm2map.tiled import resolve_tile_size

    # Tiny grid: any sane auto tile (>= 512 px) covers it both directions.
    result = resolve_tile_size(
        "auto", grid_h=10, grid_w=10, persistent_rss_bytes=10**9
    )
    assert result is None


def test_resolve_tile_size_auto_returns_int_for_large_grid():
    """``"auto"`` returns an int multiple of 256 in [512, 16384] for big outputs."""
    from isistools.csm2map.tiled import resolve_tile_size

    result = resolve_tile_size(
        "auto", grid_h=100_000, grid_w=100_000, persistent_rss_bytes=10**9
    )
    assert isinstance(result, int)
    assert 512 <= result <= 16384
    assert result % 256 == 0


def test_resolve_tile_size_invalid_raises():
    """Unknown strings and wrong types raise."""
    from isistools.csm2map.tiled import resolve_tile_size

    with pytest.raises(ValueError):
        resolve_tile_size("garbage", grid_h=10, grid_w=10, persistent_rss_bytes=10**9)
    with pytest.raises(TypeError):
        resolve_tile_size(  # type: ignore[arg-type]
            3.14, grid_h=10, grid_w=10, persistent_rss_bytes=10**9
        )


@pytest.mark.slow
def test_ctx_end_to_end():
    """End-to-end test with real CTX data. Requires ISIS + CSM deps + data."""
    pytest.skip("Requires CTX test data -- run manually")
