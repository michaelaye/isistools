"""Tiled csm2map: bound peak memory by processing the output in tile-sized windows.

For HiRISE-scale outputs (>1 Gpx) the batch path's peak memory scales
linearly with output pixel count, exceeding typical RAM. The tiled path
bounds peak at ``O(tile_size**2)`` regardless of total output size.

How it works
------------
1. Camera, shape model, and output grid are resolved exactly as in the
   batch path (see :func:`isistools.csm2map.csm2map`).
2. CSM is evaluated **once globally** on the coarse grid via
   :func:`compute_coarse_state`. This produces a tiny
   ``(nrows_coarse, ncols_coarse)`` pair of arrays (~170 KB on B17 at
   step=32) that all tiles share.
3. The output GeoTIFF is opened once with the same compression/tiling
   profile as the batch path. Tiles are processed sequentially.
   For each tile the coarse arrays are bilinear-upsampled only into
   the tile window via :func:`coordinate_map_for_window`, the input
   image is resampled into the window, and the result is written via
   ``rasterio.windows.Window``.

Tile boundaries are pixel-exact: the coarse-to-output stride in the
windowed upsample is derived from the **full** output dimensions, so
adjacent tiles agree on shared rows/cols without overlap-and-trim.

Out of scope for v1
-------------------
- CLI wiring (``csm2map`` CLI continues to use the batch path).
- Auto-dispatch threshold (caller picks the path explicitly).
- ``--clip-to-footprint`` (would re-introduce a full ``(h, w)`` bool
  via ``rasterize``, defeating the tiling).
- ``dense=True`` and ``validate=True`` flags.
- Cross-tile parallelism. Each tile uses the existing internal
  threading inside ``_resample_band`` and ``_bilinear_upsample_pair_window``.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import rasterio
from rasterio.windows import Window
from rich.console import Console

from isistools.csm2map.camera import get_image_size, load_camera
from isistools.csm2map.dem import DemRadiusSampler, resolve_shape_model
from isistools.csm2map.grid import OutputGrid, grid_from_map_file, grid_from_params
from isistools.csm2map.resample import Interpolation, resample
from isistools.csm2map.transform import (
    CoarseState,
    compute_coarse_state,
    coordinate_map_for_window,
)
from isistools.csm2map.writers import write_mapping_pvl
from isistools.io.cubes import read_isis_cube_raw

if TYPE_CHECKING:
    import numpy as np

    from isistools.csm2map.camera import TargetBody

console = Console()


# Per-pixel peak memory models for batch vs tiled paths.
#   Tile (~40 bytes/pixel): coords (8) + valid (1) + result (4) +
#     chained-&= transient (1) + scipy float64 scratch (~8) + weight
#     fragments in upsample (~4) + slack for interp-method variance.
#     Conservative — over-estimating produces smaller tiles, not OOM.
#   Batch (~25 bytes/pixel): same arrays sized to the full output;
#     no per-tile setup transients, no rasterio window overhead.
#     Calibrated against the L1 B17 measurement: peak 16.24 GB on a
#     653 M-pixel output ≈ 24.9 bytes/pixel measured, rounded up.
_BYTES_PER_TILE_PIXEL = 40
_BYTES_PER_BATCH_PIXEL = 25

# Floor / ceiling for the auto-sized tile edge.
#   - 512 px: below this, per-tile setup overhead (Python loop, np.empty,
#     thread-pool dispatch) starts to dominate the per-tile work.
#   - 16384 px: 268 M pixels @ 40 bytes ≈ 10 GB per tile, which is already
#     larger than most HiRISE-scale jobs need.
#   - Snapped to a multiple of the GeoTIFF block size (256 px) so windowed
#     writes align cleanly and we don't pay a read-modify-write penalty.
_AUTO_TILE_MIN = 512
_AUTO_TILE_MAX = 16384
_AUTO_TILE_BLOCK = 256

# OS / kernel reserve. macOS keeps ~3 GB wired + cache it can't easily
# evict; treating it as effectively-unavailable prevents the auto-sizer
# from running the system into swap. Linux is similar but smaller.
_OS_RESERVE_BYTES = 3 * 1024**3

# Headroom against `total` RAM (not `available`). On macOS, `available`
# is a conservative "free without performance hit" metric — it can
# under-report by 8-10 GB on idle 24 GB machines because compressed and
# cached pages are excluded. We use total minus persistent RSS minus
# OS reserve as the realistic work budget; 0.85 leaves a final 15%
# slack for runtime fluctuation. See the L1 measurement: idle B17
# reported `available=7.23 GB` but actually peaked at 16.24 GB without
# OOM, confirming `available` underestimates by ~10 GB.
_BUDGET_HEADROOM = 0.85


def _work_budget_bytes(persistent_rss_bytes: int) -> int:
    """How much memory the projection work can use, in bytes.

    Computed as ``(total - persistent_rss - OS reserve) * headroom``.
    Caller's persistent RSS is subtracted so the budget reflects what's
    actually available *to do new allocations*, not the OS-conservative
    ``available`` metric.
    """
    import psutil

    total = psutil.virtual_memory().total
    work = total - persistent_rss_bytes - _OS_RESERVE_BYTES
    if work <= 0:
        return 0
    return int(work * _BUDGET_HEADROOM)


def _auto_tile_size(
    *,
    grid_h: int,
    grid_w: int,
    persistent_rss_bytes: int,
    bytes_per_pixel: int = _BYTES_PER_TILE_PIXEL,
) -> int:
    """Pick a tile edge that fits the current machine's RAM budget.

    Caller measures the process's actual resident set size (post camera
    load + input read + coarse-state build). The work budget is derived
    from total RAM, the persistent baseline, and an OS reserve — see
    ``_work_budget_bytes`` for the formula. Returns a tile edge snapped
    to a 256-pixel multiple in ``[_AUTO_TILE_MIN, _AUTO_TILE_MAX]``.

    Note: callers that want batch-vs-tiled dispatch should use
    ``resolve_tile_size`` instead; that helper estimates batch peak
    first and only tiles if batch wouldn't fit.
    """
    import psutil

    total = psutil.virtual_memory().total
    available = psutil.virtual_memory().available
    budget = _work_budget_bytes(persistent_rss_bytes)
    if budget <= 0:
        return _AUTO_TILE_MIN

    max_tile_pixels = budget // bytes_per_pixel
    tile_edge = int(math.sqrt(max_tile_pixels))
    tile_edge = (tile_edge // _AUTO_TILE_BLOCK) * _AUTO_TILE_BLOCK
    tile_edge = max(_AUTO_TILE_MIN, min(tile_edge, _AUTO_TILE_MAX))

    console.print(
        f"  [dim]auto tile: {tile_edge} px "
        f"(work budget {budget / 1024**3:.2f} GB / "
        f"total {total / 1024**3:.2f} GB / "
        f"available {available / 1024**3:.2f} GB / "
        f"persistent RSS {persistent_rss_bytes / 1024**3:.2f} GB)[/dim]"
    )
    _ = grid_h, grid_w  # caller decides batch fall-through
    return tile_edge


def _batch_fits(grid_h: int, grid_w: int, persistent_rss_bytes: int) -> bool:
    """True if the batch path's estimated peak fits the work budget.

    Batch peak ≈ ``grid_h * grid_w * _BYTES_PER_BATCH_PIXEL``. If that
    fits within ``_work_budget_bytes(persistent_rss)``, the batch path
    is preferred — it's faster (no per-tile setup, no rasterio window
    overhead, single ZSTD compression pass).
    """
    budget = _work_budget_bytes(persistent_rss_bytes)
    if budget <= 0:
        return False
    estimated_batch_peak = grid_h * grid_w * _BYTES_PER_BATCH_PIXEL
    return estimated_batch_peak <= budget


def resolve_tile_size(
    spec: int | str | None,
    *,
    grid_h: int,
    grid_w: int,
    persistent_rss_bytes: int,
) -> int | None:
    """Resolve a user-supplied tile-size spec to a concrete edge or None.

    Returns ``None`` to mean "use the batch path" (single-tile or user
    forced "none"); returns an int to mean "use the tiled path with this
    edge". Callers dispatch on the return value.

    Spec semantics:
        - ``"auto"``: measure-then-decide via ``_auto_tile_size``. If the
          chosen edge covers the full output, returns ``None`` (batch).
        - ``"none"`` or ``0`` or ``None``: force batch (returns ``None``).
        - positive int: use exactly that edge. No fall-through.
    """
    if spec in (None, 0, "none", "None"):
        return None
    if isinstance(spec, str):
        if spec.lower() == "auto":
            # Prefer batch when it fits the budget — it's faster than
            # tiled (no per-tile setup, single ZSTD compression pass).
            # Only tile when batch would exceed the work budget.
            if _batch_fits(grid_h, grid_w, persistent_rss_bytes):
                budget_gb = _work_budget_bytes(persistent_rss_bytes) / 1024**3
                est_gb = grid_h * grid_w * _BYTES_PER_BATCH_PIXEL / 1024**3
                console.print(
                    f"  [dim]auto: batch path fits "
                    f"(estimated {est_gb:.2f} GB ≤ budget {budget_gb:.2f} GB)[/dim]"
                )
                return None
            edge = _auto_tile_size(
                grid_h=grid_h,
                grid_w=grid_w,
                persistent_rss_bytes=persistent_rss_bytes,
            )
            if edge >= grid_h and edge >= grid_w:
                console.print(
                    "  [dim]auto tile covers full output → "
                    "falling through to batch path[/dim]"
                )
                return None
            return edge
        try:
            spec = int(spec)
        except ValueError as e:
            msg = f"tile_size string must be 'auto', 'none', or an int (got {spec!r})"
            raise ValueError(msg) from e
    if isinstance(spec, int):
        if spec <= 0:
            return None
        return spec
    msg = f"tile_size must be int, 'auto', 'none', or None (got {type(spec).__name__})"
    raise TypeError(msg)


def project_tiled(
    *,
    state: CoarseState,
    data: "np.ndarray",
    grid: OutputGrid,
    body: "TargetBody",
    output_path: Path,
    tile_size: int,
    interpolation: Interpolation = Interpolation.BICUBIC,
    write_pvl: bool = True,
) -> Path:
    """Tile-and-write loop. Reusable core for both csm2map_tiled and ctxpipe.

    Caller is responsible for:
    - loading the camera (yielding ``state`` via compute_coarse_state)
    - reading the input image (``data``, single-band float32 ndarray)
    - building the output ``grid`` and target ``body`` description

    This function only opens the output GeoTIFF, loops tiles, resamples
    each into a window, and (optionally) writes the PVL sidecar.
    """
    if data.ndim != 2:
        msg = (
            f"project_tiled supports single-band input only "
            f"(got shape {data.shape}); use the batch path instead."
        )
        raise NotImplementedError(msg)

    n_tiles_r = (grid.height + tile_size - 1) // tile_size
    n_tiles_c = (grid.width + tile_size - 1) // tile_size
    n_tiles = n_tiles_r * n_tiles_c
    console.print(
        f"[bold]Tiling output[/bold]: {n_tiles_r} x {n_tiles_c} = "
        f"{n_tiles} tiles of up to {tile_size} x {tile_size}"
    )

    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "width": grid.width,
        "height": grid.height,
        "count": 1,
        "crs": grid.crs.to_wkt(),
        "transform": grid.transform,
        "nodata": 0.0,
        "compress": "zstd",
        "zstd_level": 3,
        "num_threads": "ALL_CPUS",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(str(output_path), "w", **profile) as dst:
        for ti, row0 in enumerate(range(0, grid.height, tile_size)):
            for tj, col0 in enumerate(range(0, grid.width, tile_size)):
                h_t = min(tile_size, grid.height - row0)
                w_t = min(tile_size, grid.width - col0)
                coord_map = coordinate_map_for_window(state, row0, col0, h_t, w_t)
                # fill_value=0.0 lets resample take its zero-transient
                # in-place fast path (result *= valid). The output GeoTIFF
                # nodata is also 0.0, so the DN value matches what the
                # batch path writes (write_geotiff does np.where(isnan)
                # to nodata=0.0 on its end).
                tile = resample(
                    data,
                    coord_map,
                    interpolation=interpolation,
                    fill_value=0.0,
                )
                dst.write(tile, 1, window=Window(col0, row0, w_t, h_t))
                console.print(
                    f"  tile {ti * n_tiles_c + tj + 1}/{n_tiles} "
                    f"@ ({row0}, {col0}) {h_t}x{w_t} written"
                )

    if write_pvl:
        pvl_path = write_mapping_pvl(output_path, grid, body)
        console.print(f"  Mapping sidecar: {pvl_path.name}")
    return output_path


def csm2map_tiled(
    input_cube: str | Path,
    output_path: str | Path,
    *,
    tile_size: int = 4096,
    # Grid definition -- either map_file OR explicit params
    map_file: str | Path | None = None,
    projection: str | None = None,
    resolution: float | None = None,
    lat_range: tuple[float, float] | None = None,
    lon_range: tuple[float, float] | None = None,
    coarse_step: int = 32,
    shape_model: str | Path | None | Literal["auto", "ellipsoid"] = "auto",
    spice_source: Literal["isis", "naif", "auto"] = "isis",
    interpolation: Interpolation = Interpolation.BICUBIC,
) -> Path:
    """Map-project an ISIS cube via the tiled output path.

    Functionally equivalent to :func:`isistools.csm2map.csm2map` on the
    pixels it writes — the per-tile coarse-grid windowed upsample is
    bit-identical to the global upsample restricted to the same window.

    Parameters
    ----------
    input_cube, output_path
        See :func:`csm2map`.
    tile_size : int
        Square tile edge in pixels. 4096 is a reasonable default; smaller
        tiles bound memory tighter but pay more per-tile setup overhead.
    map_file, projection, resolution, lat_range, lon_range
        See :func:`csm2map` for grid-definition semantics.
    coarse_step : int
        Coarse-grid spacing in pixels. Same default as :func:`csm2map`.
    shape_model, spice_source
        See :func:`csm2map`.
    interpolation : Interpolation
        Resampling method.

    Returns
    -------
    Path to the written GeoTIFF.
    """
    input_cube = Path(input_cube)
    output_path = Path(output_path)

    console.print(
        f"[bold]Loading CSM camera model[/bold] from {input_cube.name} "
        f"(SPICE source: {spice_source})"
    )
    model, body = load_camera(input_cube, spice_source=spice_source)
    n_lines, n_samples = get_image_size(model)
    console.print(f"  Input image: {n_samples} x {n_lines} (samples x lines)")
    console.print(
        f"  Target: {body.name} (NAIF {body.naif_id})  "
        f"radii eq={body.radius_equatorial_m:.1f} m  polar={body.radius_polar_m:.1f} m"
    )

    mean_radius = body.radius_mean_m
    dem_sampler: DemRadiusSampler | None = None
    if shape_model == "auto":
        dem_path = resolve_shape_model(input_cube)
        if dem_path is not None:
            dem_sampler = DemRadiusSampler(dem_path, fallback_radius=mean_radius)
            console.print(f"  Shape model: DEM {dem_path.name}")
        else:
            console.print(f"  Shape model: ellipsoid (radius {mean_radius:.1f} m)")
    elif shape_model in (None, "ellipsoid"):
        console.print(f"  Shape model: ellipsoid (radius {mean_radius:.1f} m)")
    else:
        dem_path = Path(shape_model)
        if not dem_path.exists():
            msg = f"DEM cube not found: {dem_path}"
            raise FileNotFoundError(msg)
        dem_sampler = DemRadiusSampler(dem_path, fallback_radius=mean_radius)
        console.print(f"  Shape model: DEM {dem_path.name}")

    console.print("[bold]Defining output grid[/bold]")
    grid = _build_grid(
        model=model,
        body=body,
        map_file=map_file,
        projection=projection,
        resolution=resolution,
        lat_range=lat_range,
        lon_range=lon_range,
    )
    console.print(
        f"  Output: {grid.width} x {grid.height} pixels, {grid.resolution:.2f} m/px"
    )

    console.print("[bold]Building coarse state[/bold] (one CSM pass for all tiles)")
    state = compute_coarse_state(
        model,
        grid,
        mean_radius,
        step=coarse_step,
        input_n_lines=n_lines,
        input_n_samples=n_samples,
        dem_sampler=dem_sampler,
    )
    console.print(
        f"  Coarse grid: {state.coarse_lines.shape[0]} x "
        f"{state.coarse_lines.shape[1]} CSM evaluations"
    )

    console.print(f"[bold]Reading[/bold] {input_cube.name}")
    data, _label = read_isis_cube_raw(input_cube)

    project_tiled(
        state=state,
        data=data,
        grid=grid,
        body=body,
        output_path=output_path,
        tile_size=tile_size,
        interpolation=interpolation,
    )
    console.print(f"[green bold]Done![/green bold] -> {output_path}")
    return output_path


def _build_grid(
    *,
    model,
    body,
    map_file: Path | None,
    projection: str | None,
    resolution: float | None,
    lat_range: tuple[float, float] | None,
    lon_range: tuple[float, float] | None,
) -> OutputGrid:
    """Same grid-build logic as the batch pipeline.

    Duplicated rather than refactored out of ``pipeline.py`` to keep
    the batch path's behavior provably untouched in this PR. A follow-up
    can extract a shared helper once both paths are in production.
    """
    from isistools.csm2map.camera import compute_ground_sample_distance
    from isistools.csm2map.pipeline import _derive_ground_range

    if map_file is not None:
        return grid_from_map_file(
            map_file,
            camera_lat_range=lat_range,
            camera_lon_range=lon_range,
            resolution_override=resolution,
        )

    if lat_range is None or lon_range is None:
        lat_range, lon_range = _derive_ground_range(model)

    if projection is None:
        center_lat = (lat_range[0] + lat_range[1]) / 2.0
        center_lon = (lon_range[0] + lon_range[1]) / 2.0
        projection = (
            f"+proj=eqc +lat_ts={center_lat} +lon_0={center_lon} "
            f"+a={body.radius_equatorial_m} +b={body.radius_polar_m} "
            f"+units=m +no_defs +type=crs"
        )

    if resolution is None:
        resolution = compute_ground_sample_distance(model, body)
        console.print(f"  Auto-resolution from GSD: {resolution:.4f} m/px")

    return grid_from_params(
        crs=projection,
        resolution=resolution,
        lat_min=lat_range[0],
        lat_max=lat_range[1],
        lon_min=lon_range[0],
        lon_max=lon_range[1],
    )
