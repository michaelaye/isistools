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

from pathlib import Path
from typing import Literal

import rasterio
from rasterio.windows import Window
from rich.console import Console

from isistools.csm2map.camera import get_image_size, load_camera
from isistools.csm2map.dem import DemRadiusSampler, resolve_shape_model
from isistools.csm2map.grid import OutputGrid, grid_from_map_file, grid_from_params
from isistools.csm2map.resample import Interpolation, resample
from isistools.csm2map.transform import (
    compute_coarse_state,
    coordinate_map_for_window,
)
from isistools.csm2map.writers import write_mapping_pvl
from isistools.io.cubes import read_isis_cube_raw

console = Console()


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
    if data.ndim != 2:
        msg = (
            f"Tiled path supports single-band input only "
            f"(got shape {data.shape}); use the batch csm2map() instead."
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

    pvl_path = write_mapping_pvl(output_path, grid, body)
    console.print(f"  Mapping sidecar: {pvl_path.name}")
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
