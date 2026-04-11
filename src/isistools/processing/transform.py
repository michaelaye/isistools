"""Coordinate transform: output grid pixels -> input image pixels.

Two strategies:
1. Dense: evaluate CSM groundToImage at every output pixel.  Correct but slow.
2. Coarse-grid + interpolation: evaluate on a subsampled grid, then bilinearly
   interpolate.  Much faster, with controllable accuracy.

The coarse-grid approach is essentially what ISIS ProcessRubberSheet does
(fitting bilinear patches), but implemented as a regular grid which is
trivially vectorizable and easy to reason about accuracy.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from isistools.processing.camera import ground_to_image_batch
from isistools.processing.grid import OutputGrid

if TYPE_CHECKING:
    import csmapi

    from isistools.processing.dem import DemRadiusSampler


@dataclass
class CoordinateMap:
    """Dense mapping from output pixel coords to input pixel coords."""

    input_lines: np.ndarray  # shape (height, width), float64
    input_samples: np.ndarray  # shape (height, width), float64
    valid: np.ndarray  # shape (height, width), bool

    @property
    def shape(self) -> tuple[int, int]:
        return self.input_lines.shape


def compute_transform_dense(
    model: "csmapi.RasterGM",
    grid: OutputGrid,
    surface_radius: float,
    input_n_lines: int | None = None,
    input_n_samples: int | None = None,
    dem_sampler: "DemRadiusSampler | None" = None,
) -> CoordinateMap:
    """Evaluate CSM groundToImage at every output pixel (brute-force).

    Use this for validation.  For production, use compute_transform_coarse().

    Parameters
    ----------
    model : csmapi.RasterGM
        CSM sensor model for the input image.
    grid : OutputGrid
        Fully defined output raster grid.
    surface_radius : float
        Constant fallback surface radius in meters used when no DEM is
        provided or when the DEM has nodata at a point.
    input_n_lines, input_n_samples : int, optional
        Input image dimensions for validity masking.
    dem_sampler : DemRadiusSampler, optional
        If provided, sample local body radius per point from a DEM cube
        instead of using the constant ``surface_radius``. Matches ISIS
        cam2map's behavior of using a shape model.

    Returns
    -------
    CoordinateMap
    """
    lat_rad, lon_rad = grid.pixel_to_latlon()
    if dem_sampler is not None:
        radii = dem_sampler.sample_radii(lat_rad, lon_rad)
    else:
        radii = np.full_like(lat_rad, surface_radius)

    in_lines, in_samps = ground_to_image_batch(model, lat_rad, lon_rad, radii)

    valid = np.isfinite(in_lines) & np.isfinite(in_samps)

    if input_n_lines is not None and input_n_samples is not None:
        in_image = (
            (in_lines >= -0.5)
            & (in_lines <= input_n_lines - 0.5)
            & (in_samps >= -0.5)
            & (in_samps <= input_n_samples - 0.5)
        )
        valid &= in_image

    return CoordinateMap(
        input_lines=in_lines,
        input_samples=in_samps,
        valid=valid,
    )


def compute_transform_coarse(
    model: "csmapi.RasterGM",
    grid: OutputGrid,
    surface_radius: float,
    step: int = 16,
    input_n_lines: int | None = None,
    input_n_samples: int | None = None,
    dem_sampler: "DemRadiusSampler | None" = None,
) -> CoordinateMap:
    """Evaluate CSM on a coarse grid, then bilinearly interpolate.

    This is the production path.  The coarse grid is evaluated at every
    ``step``-th pixel, then scipy RegularGridInterpolator fills in the rest.

    Parameters
    ----------
    model : csmapi.RasterGM
        CSM sensor model.
    grid : OutputGrid
        Output raster grid definition.
    surface_radius : float
        Constant fallback surface radius in meters used when no DEM is
        provided or when the DEM has nodata at a point.
    step : int
        Coarse grid spacing in pixels.  Smaller = more accurate, slower.
        16 is a good default for CTX; use 8 for HiRISE with high-res DEM.
    input_n_lines, input_n_samples : int, optional
        Input image dimensions for validity masking.  If provided, points
        that fall outside the input image are masked as invalid.
    dem_sampler : DemRadiusSampler, optional
        If provided, sample local body radius per coarse-grid point from
        a DEM cube instead of using the constant ``surface_radius``.

    Returns
    -------
    CoordinateMap
    """
    from pyproj import CRS, Transformer

    h, w = grid.height, grid.width

    # Build coarse grid indices
    rows_coarse = np.arange(0, h, step)
    cols_coarse = np.arange(0, w, step)
    # Ensure we include the last row/col
    if rows_coarse[-1] != h - 1:
        rows_coarse = np.append(rows_coarse, h - 1)
    if cols_coarse[-1] != w - 1:
        cols_coarse = np.append(cols_coarse, w - 1)

    cc, rr = np.meshgrid(cols_coarse, rows_coarse)

    # Pixel center -> map coords via affine
    x = grid.transform.c + (cc + 0.5) * grid.transform.a
    y = grid.transform.f + (rr + 0.5) * grid.transform.e

    # Map coords -> lat/lon via inverse projection
    transformer = Transformer.from_crs(grid.crs, CRS.from_epsg(4326), always_xy=True)
    lon_deg, lat_deg = transformer.transform(x, y)
    lat_rad = np.deg2rad(lat_deg)
    lon_rad = np.deg2rad(lon_deg)
    if dem_sampler is not None:
        radii = dem_sampler.sample_radii(lat_rad, lon_rad)
    else:
        radii = np.full_like(lat_rad, surface_radius)

    # Evaluate CSM on coarse grid
    coarse_lines, coarse_samps = ground_to_image_batch(model, lat_rad, lon_rad, radii)

    # Do NOT mask out-of-bounds coarse points. Masking them causes
    # RegularGridInterpolator to propagate NaN to the entire 16x16 cell
    # via bilinear interpolation, which produces a ragged, ~1-pixel-
    # inaccurate edge. Instead, let the mapping extrapolate smoothly and
    # apply the per-pixel bounds check AFTER interpolation.

    # Interpolate to full resolution
    interp_lines = RegularGridInterpolator(
        (rows_coarse.astype(float), cols_coarse.astype(float)),
        coarse_lines,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    interp_samps = RegularGridInterpolator(
        (rows_coarse.astype(float), cols_coarse.astype(float)),
        coarse_samps,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )

    # Evaluate at all output pixels
    all_rows = np.arange(h)
    all_cols = np.arange(w)
    cc_full, rr_full = np.meshgrid(all_cols, all_rows)
    pts = np.column_stack([rr_full.ravel(), cc_full.ravel()])

    in_lines = interp_lines(pts).reshape(h, w)
    in_samps = interp_samps(pts).reshape(h, w)

    valid = np.isfinite(in_lines) & np.isfinite(in_samps)

    # Per-pixel bounds check against the input image extent.
    if input_n_lines is not None and input_n_samples is not None:
        in_image = (
            (in_lines >= -0.5)
            & (in_lines <= input_n_lines - 0.5)
            & (in_samps >= -0.5)
            & (in_samps <= input_n_samples - 0.5)
        )
        valid &= in_image

    return CoordinateMap(
        input_lines=in_lines,
        input_samples=in_samps,
        valid=valid,
    )


def validate_coarse_vs_dense(
    model: "csmapi.RasterGM",
    grid: OutputGrid,
    coord_map: CoordinateMap,
    surface_radius: float,
    n_check: int = 1000,
    tolerance: float = 0.5,
) -> dict:
    """Spot-check coarse-interpolated transform against dense CSM evaluation.

    Randomly samples ``n_check`` output pixels, evaluates CSM directly, and
    compares against the interpolated coordinate map.

    Returns
    -------
    dict with keys:
        max_error_line, max_error_sample : float
        mean_error_line, mean_error_sample : float
        n_checked : int
        n_failed : int  (error > tolerance)
    """
    from pyproj import CRS, Transformer

    h, w = coord_map.shape
    rng = np.random.default_rng(42)

    # Sample random valid pixels
    valid_idx = np.argwhere(coord_map.valid)
    if len(valid_idx) < n_check:
        n_check = len(valid_idx)
    chosen = rng.choice(len(valid_idx), size=n_check, replace=False)
    check_rows = valid_idx[chosen, 0]
    check_cols = valid_idx[chosen, 1]

    # Compute lat/lon for these pixels
    x = grid.transform.c + (check_cols + 0.5) * grid.transform.a
    y = grid.transform.f + (check_rows + 0.5) * grid.transform.e

    transformer = Transformer.from_crs(grid.crs, CRS.from_epsg(4326), always_xy=True)
    lon_deg, lat_deg = transformer.transform(x, y)

    radii = np.full_like(lat_deg, surface_radius)
    true_lines, true_samps = ground_to_image_batch(
        model, np.deg2rad(lat_deg), np.deg2rad(lon_deg), radii
    )

    interp_lines = coord_map.input_lines[check_rows, check_cols]
    interp_samps = coord_map.input_samples[check_rows, check_cols]

    err_l = np.abs(true_lines - interp_lines)
    err_s = np.abs(true_samps - interp_samps)

    valid_check = np.isfinite(err_l) & np.isfinite(err_s)
    err_l = err_l[valid_check]
    err_s = err_s[valid_check]

    return {
        "max_error_line": float(np.max(err_l)) if len(err_l) > 0 else 0.0,
        "max_error_sample": float(np.max(err_s)) if len(err_s) > 0 else 0.0,
        "mean_error_line": float(np.mean(err_l)) if len(err_l) > 0 else 0.0,
        "mean_error_sample": float(np.mean(err_s)) if len(err_s) > 0 else 0.0,
        "n_checked": int(np.sum(valid_check)),
        "n_failed": int(np.sum((err_l > tolerance) | (err_s > tolerance))),
    }
