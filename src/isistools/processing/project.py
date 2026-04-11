"""Main cam2map pipeline: orchestrate camera, grid, transform, resample, I/O."""

from pathlib import Path
from typing import Literal

import numpy as np
from rich.console import Console

from isistools.io.cubes import read_isis_cube_raw
from isistools.processing.camera import get_image_size, get_target_radii, load_camera
from isistools.processing.dem import DemRadiusSampler, resolve_shape_model
from isistools.processing.grid import OutputGrid, grid_from_map_file, grid_from_params
from isistools.processing.resample import Interpolation, resample
from isistools.processing.transform import (
    compute_transform_coarse,
    compute_transform_dense,
    validate_coarse_vs_dense,
)
from isistools.processing.writers import write_geotiff

console = Console()


def project(
    input_cube: str | Path,
    output_path: str | Path,
    *,
    # Grid definition -- either map_file OR explicit params
    map_file: str | Path | None = None,
    projection: str | None = None,
    resolution: float | None = None,
    lat_range: tuple[float, float] | None = None,
    lon_range: tuple[float, float] | None = None,
    # Transform options
    coarse_step: int = 16,
    dense: bool = False,
    validate: bool = False,
    clip_to_footprint: bool = False,
    shape_model: str | Path | None | Literal["auto", "ellipsoid"] = "auto",
    # Resampling
    interpolation: Interpolation = Interpolation.BICUBIC,
    # Output
    output_format: Literal["geotiff", "cube"] = "geotiff",
) -> Path:
    """Map-project an ISIS cube using CSM camera model.

    Parameters
    ----------
    input_cube : path-like
        Path to spiceinit'd Level 1 ISIS cube.
    output_path : path-like
        Output file path.
    map_file : path-like, optional
        ISIS MAP file for grid definition (cam2map-compatible).
    projection : str, optional
        PROJ string or projection name (if not using map_file).
    resolution : float, optional
        Pixel resolution in meters/pixel.
    lat_range, lon_range : tuple, optional
        (min, max) ground range in degrees.
    coarse_step : int
        Coarse grid spacing for the interpolation approach.
    dense : bool
        If True, evaluate CSM at every pixel (slow, for validation).
    validate : bool
        If True, spot-check the coarse transform against dense evaluation.
    clip_to_footprint : bool
        If True, clip output to the footprint polygon stored in the cube
        by ``footprintinit``. This reproduces ISIS cam2map's behavior but
        inherits ``footprintinit``'s precision limits. Default False -
        our camera-model-based mask is more accurate than the polygon.
        Use this flag only when comparing against ISIS cam2map output.
    shape_model : str, Path, "auto", "ellipsoid", or None
        Shape model used for the body's local radius:
          - ``"auto"`` (default): read ``Kernels.ShapeModel`` from the
            input cube label and use that DEM if present, otherwise
            fall back to ellipsoid. Matches ISIS cam2map's default.
          - ``"ellipsoid"`` or ``None``: use the constant mean radius
            from the cube's target body.
          - ``Path``: explicit path to an ISIS DEM cube.
    interpolation : Interpolation
        Pixel interpolation method.
    output_format : str
        Output format ("geotiff" or "cube").

    Returns
    -------
    Path to the output file.
    """
    input_cube = Path(input_cube)
    output_path = Path(output_path)

    # Step 1: Load camera model
    console.print(f"[bold]Loading CSM camera model[/bold] from {input_cube.name}")
    model = load_camera(input_cube)
    n_lines, n_samples = get_image_size(model)
    console.print(f"  Input image: {n_samples} x {n_lines} (samples x lines)")

    # Step 2: Get target radii for sphere approximation
    eq_r, polar_r = get_target_radii(input_cube)
    # Mean radius used as fallback when no DEM is given or DEM has nodata
    mean_radius = (2 * eq_r + polar_r) / 3.0

    # Step 2b: Resolve shape model
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

    # Step 3: Define output grid
    console.print("[bold]Defining output grid[/bold]")
    grid = _build_grid(
        model=model,
        input_cube=input_cube,
        map_file=map_file,
        projection=projection,
        resolution=resolution,
        lat_range=lat_range,
        lon_range=lon_range,
    )
    console.print(f"  Output: {grid.width} x {grid.height} pixels, {grid.resolution:.2f} m/px")
    console.print(
        f"  Lat: [{grid.lat_min:.4f}, {grid.lat_max:.4f}]  "
        f"Lon: [{grid.lon_min:.4f}, {grid.lon_max:.4f}]"
    )

    # Step 4: Compute coordinate transform
    if dense:
        console.print("[bold]Computing dense transform[/bold] (every pixel)...")
        coord_map = compute_transform_dense(
            model,
            grid,
            mean_radius,
            input_n_lines=n_lines,
            input_n_samples=n_samples,
            dem_sampler=dem_sampler,
        )
    else:
        n_coarse = (grid.height // coarse_step + 1) * (grid.width // coarse_step + 1)
        console.print(
            f"[bold]Computing coarse transform[/bold] "
            f"(step={coarse_step}, ~{n_coarse:,} CSM calls)..."
        )
        coord_map = compute_transform_coarse(
            model,
            grid,
            mean_radius,
            step=coarse_step,
            input_n_lines=n_lines,
            input_n_samples=n_samples,
            dem_sampler=dem_sampler,
        )

    # Step 4b: Optional footprint-polygon clipping (ISIS cam2map compat mode)
    if clip_to_footprint:
        console.print("[bold]Clipping to footprint polygon[/bold] (ISIS compat mode)")
        polygon_mask = _rasterize_footprint(input_cube, grid)
        coord_map.valid &= polygon_mask

    n_valid = int(np.sum(coord_map.valid))
    n_total = grid.height * grid.width
    console.print(f"  Valid pixels: {n_valid:,} / {n_total:,} ({100 * n_valid / n_total:.1f}%)")

    # Step 4c: Optional validation
    if validate and not dense:
        console.print("[bold]Validating[/bold] coarse transform (1000 random points)...")
        stats = validate_coarse_vs_dense(model, grid, coord_map, mean_radius)
        console.print(
            f"  Max error: {stats['max_error_line']:.4f} lines, "
            f"{stats['max_error_sample']:.4f} samples"
        )
        console.print(
            f"  Mean error: {stats['mean_error_line']:.4f} lines, "
            f"{stats['mean_error_sample']:.4f} samples"
        )
        if stats["n_failed"] > 0:
            console.print(
                f"  [yellow]Warning: {stats['n_failed']} / {stats['n_checked']} "
                f"points exceeded 0.5 pixel tolerance[/yellow]"
            )

    # Step 5: Read input image
    console.print(f"[bold]Reading[/bold] {input_cube.name}")
    data, _label = read_isis_cube_raw(input_cube)

    # Step 6: Resample
    console.print(f"[bold]Resampling[/bold] ({interpolation.value})")
    projected = resample(data, coord_map, interpolation=interpolation, fill_value=np.nan)

    # Step 7: Write output
    console.print(f"[bold]Writing[/bold] {output_path.name}")
    if output_format == "geotiff":
        result = write_geotiff(output_path, projected, grid, nodata=0.0)
    else:
        msg = "ISIS cube output not yet implemented"
        raise NotImplementedError(msg)

    console.print(f"[green bold]Done![/green bold] -> {result}")
    return result


def _rasterize_footprint(input_cube: Path, grid: OutputGrid) -> np.ndarray:
    """Rasterize the cube's footprint polygon onto the output grid.

    Reads the footprint polygon stored in the cube (by ``footprintinit``)
    and returns a boolean mask with True inside the polygon. This matches
    ISIS cam2map's behavior of clipping output to the stored footprint.

    Parameters
    ----------
    input_cube : Path
        Path to the ISIS cube containing a Polygon blob.
    grid : OutputGrid
        Output grid defining the mask shape and transform.

    Returns
    -------
    ndarray of bool, shape (height, width)
    """
    from pyproj import CRS, Transformer
    from rasterio.features import rasterize
    from shapely.ops import transform as shapely_transform

    from isistools.io.footprints import read_footprint

    geom = read_footprint(input_cube)

    # Project polygon from lon/lat (EPSG:4326) to the output CRS
    tr = Transformer.from_crs(CRS.from_epsg(4326), grid.crs, always_xy=True)
    poly_in_map = shapely_transform(lambda x, y: tr.transform(x, y), geom)

    mask = rasterize(
        [(poly_in_map, 1)],
        out_shape=(grid.height, grid.width),
        transform=grid.transform,
        fill=0,
        dtype=np.uint8,
    )
    return mask.astype(bool)


def _build_grid(
    model,
    input_cube: Path,
    map_file: Path | None,
    projection: str | None,
    resolution: float | None,
    lat_range: tuple[float, float] | None,
    lon_range: tuple[float, float] | None,
) -> OutputGrid:
    """Build output grid from MAP file or explicit parameters."""
    if map_file is not None:
        return grid_from_map_file(
            map_file,
            camera_lat_range=lat_range,
            camera_lon_range=lon_range,
            resolution_override=resolution,
        )

    if projection is None:
        # Default to equirectangular with Mars radii
        projection = (
            "+proj=eqc +lat_ts=0 +lon_0=0 +a=3396190 +b=3376200 +units=m +no_defs +type=crs"
        )

    if lat_range is None or lon_range is None:
        # Try to derive from the camera model image corners
        lat_range, lon_range = _derive_ground_range(model)

    if resolution is None:
        msg = "Must specify resolution (meters/pixel) or use a MAP file"
        raise ValueError(msg)

    return grid_from_params(
        crs=projection,
        resolution=resolution,
        lat_min=lat_range[0],
        lat_max=lat_range[1],
        lon_min=lon_range[0],
        lon_max=lon_range[1],
    )


def _derive_ground_range(model) -> tuple[tuple[float, float], tuple[float, float]]:
    """Derive lat/lon ground range by probing CSM model at image corners + edges.

    Returns (lat_min, lat_max), (lon_min, lon_max) in degrees.
    """
    import csmapi

    size = model.getImageSize()
    n_lines = size.line
    n_samps = size.samp

    # Sample corners and edge midpoints
    probes = [
        (0.5, 0.5),
        (0.5, n_samps - 0.5),
        (n_lines - 0.5, 0.5),
        (n_lines - 0.5, n_samps - 0.5),
        (n_lines / 2, 0.5),
        (n_lines / 2, n_samps - 0.5),
        (0.5, n_samps / 2),
        (n_lines - 0.5, n_samps / 2),
        (n_lines / 2, n_samps / 2),
    ]

    # Also sample along the edges at ~100 points
    for frac in np.linspace(0, 1, 100):
        line = 0.5 + frac * (n_lines - 1)
        probes.append((line, 0.5))
        probes.append((line, n_samps - 0.5))
    for frac in np.linspace(0, 1, 50):
        samp = 0.5 + frac * (n_samps - 1)
        probes.append((0.5, samp))
        probes.append((n_lines - 0.5, samp))

    lats = []
    lons = []
    for line, samp in probes:
        try:
            ip = csmapi.ImageCoord(line, samp)
            gp = model.imageToGround(ip, 0.0)  # height=0
            # ECEF -> lat/lon
            r = np.sqrt(gp.x**2 + gp.y**2 + gp.z**2)
            lat = np.degrees(np.arcsin(gp.z / r))
            lon = np.degrees(np.arctan2(gp.y, gp.x))
            lats.append(lat)
            lons.append(lon)
        except Exception:
            continue

    if not lats:
        msg = "Could not derive ground range from camera model"
        raise RuntimeError(msg)

    # Add a small buffer (1%)
    lat_span = max(lats) - min(lats)
    lon_span = max(lons) - min(lons)
    buf_lat = lat_span * 0.01
    buf_lon = lon_span * 0.01

    return (
        (min(lats) - buf_lat, max(lats) + buf_lat),
        (min(lons) - buf_lon, max(lons) + buf_lon),
    )
