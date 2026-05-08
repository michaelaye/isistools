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

from isistools.csm2map.camera import ground_to_image_batch
from isistools.csm2map.grid import OutputGrid

if TYPE_CHECKING:
    import csmapi

    from isistools.csm2map.dem import DemRadiusSampler


def _bilinear_upsample_pair(
    coarse_a: np.ndarray,
    coarse_b: np.ndarray,
    h: int,
    w: int,
    workers: int | None = None,
) -> np.ndarray:
    """Bilinearly upsample two co-aligned coarse 2D arrays onto (h, w).

    Returns a single ``(2, h, w)`` float32 array. ``out[0]`` is the
    upsample of ``coarse_a``, ``out[1]`` is ``coarse_b``. Caller can
    pass this straight into ``CoordinateMap.coords`` — downstream slices
    (e.g. per-stripe in resample) are then plain views, not copies.

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

    # Output as a single (2, h, w) buffer. out[0] / out[1] are views;
    # downstream slicing (per-stripe in resample) avoids any np.stack.
    out = np.empty((2, h, w), dtype=np.float32)
    out_a = out[0]
    out_b = out[1]

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

    return out


@dataclass
class CoordinateMap:
    """Dense mapping from output pixel coords to input pixel coords.

    ``coords[0]`` holds input-line values, ``coords[1]`` holds
    input-sample values, both ``(h, w)`` float32. Storing them as one
    contiguous ``(2, h, w)`` buffer lets per-stripe consumers (the
    resampler) slice with ``coords[:, r0:r1]`` — a view, no copy —
    instead of stacking them per worker. ``input_lines`` and
    ``input_samples`` remain available as zero-cost property views for
    callers that only need one channel.
    """

    coords: np.ndarray  # shape (2, height, width)
    valid: np.ndarray  # shape (height, width), bool

    @property
    def input_lines(self) -> np.ndarray:
        return self.coords[0]

    @property
    def input_samples(self) -> np.ndarray:
        return self.coords[1]

    @property
    def shape(self) -> tuple[int, int]:
        return self.coords.shape[1:]


@dataclass
class CoarseState:
    """Coarse-grid CSM evaluation, sufficient to upsample any output window.

    For the tiled output path: evaluate CSM ONCE on the coarse grid that
    spans the full output, then upsample only into per-tile windows.
    The coarse arrays are tiny (~170 KB for a CTX 12k × 53k output at
    step=32), so holding them resident across all tiles is free.

    Attributes
    ----------
    coarse_lines, coarse_samps
        Output of ``ground_to_image_batch`` on the coarse grid; shape
        ``(nrows_coarse, ncols_coarse)``, dtype matching what
        ``ground_to_image_batch`` returns (currently float32).
    grid_h, grid_w
        Full output dimensions. Used to derive coarse-to-output stride
        so that any window upsamples consistently with a single global
        upsample of the same coarse arrays.
    input_n_lines, input_n_samples
        Input image extent for the in-bounds validity check; mirrors
        ``compute_transform_coarse``'s arguments.
    """

    coarse_lines: np.ndarray
    coarse_samps: np.ndarray
    grid_h: int
    grid_w: int
    input_n_lines: int | None
    input_n_samples: int | None


def _bilinear_upsample_pair_window(
    coarse_a: np.ndarray,
    coarse_b: np.ndarray,
    full_h: int,
    full_w: int,
    row0: int,
    col0: int,
    h_t: int,
    w_t: int,
    workers: int | None = None,
) -> np.ndarray:
    """Bilinear-upsample a (h_t, w_t) window starting at (row0, col0).

    The coarse-to-output stride is derived from ``full_h`` / ``full_w``
    (NOT from the window size), so a window upsample produces values
    bit-identical to slicing the same window out of a global upsample
    on ``(full_h, full_w)``. Tile boundaries are pixel-exact — no
    overlap, no seam-trim required.

    Returns a single ``(2, h_t, w_t)`` float32 buffer suitable for
    ``CoordinateMap.coords``.
    """
    nrc, ncc = coarse_a.shape
    step_r = (full_h - 1) / (nrc - 1) if nrc > 1 else 1.0
    step_c = (full_w - 1) / (ncc - 1) if ncc > 1 else 1.0

    a32 = coarse_a.astype(np.float32, copy=False)
    b32 = coarse_b.astype(np.float32, copy=False)

    # Window-local column indices in coarse-grid coordinates
    c = np.arange(col0, col0 + w_t, dtype=np.float32) / np.float32(step_c)
    c0 = np.floor(c).astype(np.intp)
    np.clip(c0, 0, ncc - 2, out=c0)
    c1 = c0 + 1
    cf = (c - c0.astype(np.float32))[None, :]
    one_cf = np.float32(1.0) - cf
    c0_row = c0[None, :]
    c1_row = c1[None, :]

    out = np.empty((2, h_t, w_t), dtype=np.float32)
    out_a = out[0]
    out_b = out[1]

    def _process_stripe(r_lo: int, r_hi: int) -> None:
        # r_lo/r_hi are window-local; convert to global rows for stride
        r = np.arange(row0 + r_lo, row0 + r_hi, dtype=np.float32) / np.float32(step_r)
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

    if workers <= 1 or h_t < 256:
        _process_stripe(0, h_t)
    else:
        stripe = (h_t + workers - 1) // workers
        ranges = [(i * stripe, min((i + 1) * stripe, h_t)) for i in range(workers)]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda r: _process_stripe(*r), ranges))

    return out


def compute_coarse_state(
    model: "csmapi.RasterGM",
    grid: OutputGrid,
    surface_radius: float,
    step: int = 16,
    input_n_lines: int | None = None,
    input_n_samples: int | None = None,
    dem_sampler: "DemRadiusSampler | None" = None,
) -> CoarseState:
    """Evaluate CSM on the coarse grid for a full output, no upsample.

    This is the global setup step for the tiled path. Each subsequent
    tile calls ``coordinate_map_for_window`` with the same
    ``CoarseState``; CSM is therefore evaluated exactly once across all
    tiles.
    """
    from pyproj import CRS, Transformer

    h, w = grid.height, grid.width

    nrows_coarse = max(2, (h - 1) // step + 2)
    ncols_coarse = max(2, (w - 1) // step + 2)
    rows_coarse = np.linspace(0.0, h - 1, nrows_coarse)
    cols_coarse = np.linspace(0.0, w - 1, ncols_coarse)

    cc, rr = np.meshgrid(cols_coarse, rows_coarse)

    x = grid.transform.c + (cc + 0.5) * grid.transform.a
    y = grid.transform.f + (rr + 0.5) * grid.transform.e

    transformer = Transformer.from_crs(grid.crs, CRS.from_epsg(4326), always_xy=True)
    lon_deg, lat_deg = transformer.transform(x, y)
    lat_rad = np.deg2rad(lat_deg)
    lon_rad = np.deg2rad(lon_deg)
    if dem_sampler is not None:
        radii = dem_sampler.sample_radii(lat_rad, lon_rad)
    else:
        radii = np.full_like(lat_rad, surface_radius)

    coarse_lines, coarse_samps = ground_to_image_batch(model, lat_rad, lon_rad, radii)

    return CoarseState(
        coarse_lines=coarse_lines,
        coarse_samps=coarse_samps,
        grid_h=h,
        grid_w=w,
        input_n_lines=input_n_lines,
        input_n_samples=input_n_samples,
    )


def coordinate_map_for_window(
    state: CoarseState,
    row0: int,
    col0: int,
    h_t: int,
    w_t: int,
) -> CoordinateMap:
    """Build a CoordinateMap for one output tile from a global CoarseState.

    Bilinear-upsamples only into the (h_t, w_t) window starting at
    (row0, col0) in the full output, then constructs the validity mask
    via the same chained &= logic as ``compute_transform_coarse``.
    """
    coords = _bilinear_upsample_pair_window(
        state.coarse_lines,
        state.coarse_samps,
        state.grid_h,
        state.grid_w,
        row0,
        col0,
        h_t,
        w_t,
    )
    in_lines = coords[0]
    in_samps = coords[1]

    valid = np.isfinite(in_lines)
    valid &= np.isfinite(in_samps)

    if state.input_n_lines is not None and state.input_n_samples is not None:
        valid &= in_lines >= -0.5
        valid &= in_lines <= state.input_n_lines - 0.5
        valid &= in_samps >= -0.5
        valid &= in_samps <= state.input_n_samples - 0.5

    return CoordinateMap(coords=coords, valid=valid)


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

    # Build validity in place via chained &= so the (h,w) bool
    # intermediates from the bounds check don't co-exist.
    valid = np.isfinite(in_lines)
    valid &= np.isfinite(in_samps)

    if input_n_lines is not None and input_n_samples is not None:
        valid &= in_lines >= -0.5
        valid &= in_lines <= input_n_lines - 0.5
        valid &= in_samps >= -0.5
        valid &= in_samps <= input_n_samples - 0.5

    coords = np.empty((2, *in_lines.shape), dtype=in_lines.dtype)
    coords[0] = in_lines
    coords[1] = in_samps
    return CoordinateMap(coords=coords, valid=valid)


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
    # via a single shared-weight pass. Returns a single (2, h, w) buffer
    # so per-stripe consumers can slice with no copy. See
    # _bilinear_upsample_pair for the memory-bandwidth reasoning.
    coords = _bilinear_upsample_pair(coarse_lines, coarse_samps, h, w)
    in_lines = coords[0]
    in_samps = coords[1]

    # Build validity in place via chained &= so the (h,w) bool
    # intermediates from the bounds check don't co-exist.
    valid = np.isfinite(in_lines)
    valid &= np.isfinite(in_samps)

    # Per-pixel bounds check against the input image extent.
    if input_n_lines is not None and input_n_samples is not None:
        valid &= in_lines >= -0.5
        valid &= in_lines <= input_n_lines - 0.5
        valid &= in_samps >= -0.5
        valid &= in_samps <= input_n_samples - 0.5

    return CoordinateMap(coords=coords, valid=valid)


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
