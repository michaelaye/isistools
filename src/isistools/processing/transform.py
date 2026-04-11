"""Coordinate transform: output grid pixels -> input image pixels.

Two strategies:
1. Dense: evaluate CSM groundToImage at every output pixel.  Correct but slow.
2. Coarse-grid + interpolation: evaluate on a subsampled grid, then bilinearly
   interpolate.  Much faster, with controllable accuracy.

The coarse-grid approach is essentially what ISIS ProcessRubberSheet does
(fitting bilinear patches), but implemented as a regular grid which is
trivially vectorizable and easy to reason about accuracy.
"""

import os
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from isistools.processing.camera import ground_to_image_batch
from isistools.processing.grid import OutputGrid

if TYPE_CHECKING:
    import csmapi

    from isistools.processing.dem import DemRadiusSampler


def _bilinear_upsample_pair(
    coarse_a: np.ndarray,
    coarse_b: np.ndarray,
    h: int,
    w: int,
    workers: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Bilinearly upsample two co-aligned coarse 2D arrays onto (h, w).

    ``coarse_a`` and ``coarse_b`` must have the same shape and both
    describe a strictly uniform coarse grid spanning pixel (0, 0) to
    pixel (h-1, w-1). This is enforced by ``compute_transform_coarse``
    via ``np.linspace``.

    Processing both arrays together is critical for speed: the
    bottleneck is memory bandwidth on the 46M-pixel output, and the
    per-row/per-col index/fraction arrays are shared across the two
    outputs, so we amortize the index computation.

    All arithmetic is done in ``float32`` (halves the memory traffic
    vs ``float64``). The coarse line/sample values easily fit — CTX
    is at most 12,288 lines, which float32 represents exactly.

    The computation is threaded across horizontal stripes of the
    output. Numpy broadcast operations on large float32 arrays
    release the GIL, so this gives near-linear speedup up to the
    memory-bandwidth limit (~4x at 4 workers on an 8-core machine).
    """
    nrc, ncc = coarse_a.shape
    step_r = (h - 1) / (nrc - 1) if nrc > 1 else 1.0
    step_c = (w - 1) / (ncc - 1) if ncc > 1 else 1.0

    # Downcast the coarse arrays to float32 once (they're tiny, ~170 KB
    # each, so this is free).
    a32 = coarse_a.astype(np.float32, copy=False)
    b32 = coarse_b.astype(np.float32, copy=False)

    # Continuous coarse-grid column indices for every output column
    c = np.arange(w, dtype=np.float32) / np.float32(step_c)
    c0 = np.floor(c).astype(np.intp)
    np.clip(c0, 0, ncc - 2, out=c0)
    c1 = c0 + 1
    cf = (c - c0.astype(np.float32))[None, :]
    one_cf = np.float32(1.0) - cf
    c0_row = c0[None, :]
    c1_row = c1[None, :]

    # Output arrays — allocated once, filled per stripe in workers
    out_a = np.empty((h, w), dtype=np.float32)
    out_b = np.empty((h, w), dtype=np.float32)

    def _process_stripe(r_lo: int, r_hi: int) -> None:
        r = np.arange(r_lo, r_hi, dtype=np.float32) / np.float32(step_r)
        r0 = np.floor(r).astype(np.intp)
        np.clip(r0, 0, nrc - 2, out=r0)
        r1 = r0 + 1

        rf = (r - r0.astype(np.float32))[:, None]
        one_rf = np.float32(1.0) - rf

        w00 = one_rf * one_cf
        w01 = one_rf * cf
        w10 = rf * one_cf
        w11 = rf * cf

        r0c = r0[:, None]
        r1c = r1[:, None]

        # Gather corners for channel A and combine
        out_a[r_lo:r_hi] = (
            a32[r0c, c0_row] * w00
            + a32[r0c, c1_row] * w01
            + a32[r1c, c0_row] * w10
            + a32[r1c, c1_row] * w11
        )
        out_b[r_lo:r_hi] = (
            b32[r0c, c0_row] * w00
            + b32[r0c, c1_row] * w01
            + b32[r1c, c0_row] * w10
            + b32[r1c, c1_row] * w11
        )

    if workers is None:
        workers = max(1, (os.cpu_count() or 1) // 2)

    if workers <= 1 or h < 256:
        _process_stripe(0, h)
    else:
        stripe = (h + workers - 1) // workers
        ranges = [(i * stripe, min((i + 1) * stripe, h)) for i in range(workers)]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda r: _process_stripe(*r), ranges))

    return out_a, out_b


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

    # Build a strictly uniform coarse grid that starts at 0 and ends at
    # (h - 1, w - 1). Using linspace (instead of arange + manually
    # appending h-1) guarantees even spacing, which lets us replace the
    # slow RegularGridInterpolator with a vectorized numpy bilinear
    # interpolation downstream.
    nrows_coarse = max(2, (h - 1) // step + 2)
    ncols_coarse = max(2, (w - 1) // step + 2)
    rows_coarse = np.linspace(0.0, h - 1, nrows_coarse)
    cols_coarse = np.linspace(0.0, w - 1, ncols_coarse)

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
    # interpolation to propagate NaN to the entire coarse cell, which
    # produces a ragged, ~1-pixel-inaccurate edge. Instead, let the
    # mapping extrapolate smoothly and apply the per-pixel bounds check
    # AFTER interpolation.

    # Upsample the coarse line/sample arrays to the full (h, w) grid
    # via a single shared-weight pass. See _bilinear_upsample_pair for
    # the memory-bandwidth reasoning.
    in_lines, in_samps = _bilinear_upsample_pair(coarse_lines, coarse_samps, h, w)

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
