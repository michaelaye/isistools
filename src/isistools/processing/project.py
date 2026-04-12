"""Main csm2map pipeline: orchestrate camera, grid, transform, resample, I/O."""

import time
from contextlib import contextmanager
from pathlib import Path
from typing import Literal

import numpy as np
from rich.console import Console

from isistools.io.cubes import read_isis_cube_raw
from isistools.processing.camera import (
    TargetBody,
    compute_ground_sample_distance,
    get_image_size,
    load_camera,
)
from isistools.processing.dem import DemRadiusSampler, resolve_shape_model
from isistools.processing.grid import OutputGrid, grid_from_map_file, grid_from_params
from isistools.processing.resample import Interpolation, resample
from isistools.processing.transform import (
    compute_transform_coarse,
    compute_transform_dense,
    validate_coarse_vs_dense,
)
from isistools.processing.writers import write_geotiff, write_mapping_pvl

console = Console()


@contextmanager
def _stage(timings: dict | None, name: str):
    """Record wall time of a pipeline stage into ``timings`` if not None."""
    if timings is None:
        yield
        return
    t0 = time.perf_counter()
    try:
        yield
    finally:
        timings[name] = time.perf_counter() - t0


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
    coarse_step: int = 32,
    dense: bool = False,
    validate: bool = False,
    clip_to_footprint: bool = False,
    shape_model: str | Path | None | Literal["auto", "ellipsoid"] = "auto",
    spice_source: Literal["isis", "naif", "auto"] = "isis",
    # Resampling
    interpolation: Interpolation = Interpolation.BICUBIC,
    # Output
    output_format: Literal["geotiff", "cube"] = "geotiff",
    # Profiling
    profile: bool = False,
) -> Path:
    """Map-project an ISIS cube using CSM camera model.

    Parameters
    ----------
    input_cube : path-like
        Path to spiceinit'd Level 1 ISIS cube.
    output_path : path-like
        Output file path.
    map_file : path-like, optional
        ISIS MAP file for grid definition (ISIS cam2map-compatible).
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
        If True, apply an extra mask from the footprint polygon stored in
        the cube by ``footprintinit``. This does NOT reproduce ISIS
        ``cam2map`` behavior — empirical testing (see
        ``docs/csm2map-design.md``) showed that ``cam2map`` ignores the
        polygon entirely; its output is determined purely by the camera
        model and shape model. The polygon mask is therefore strictly
        *additional* to our camera-model mask and can only remove valid
        pixels that were otherwise correctly projected. Kept as an
        escape hatch for callers who want a polygon-clipped output for
        their own downstream reasons; not recommended for comparison
        against ``cam2map``. Default False.
    shape_model : str, Path, "auto", "ellipsoid", or None
        Shape model used for the body's local radius:
          - ``"auto"`` (default): read ``Kernels.ShapeModel`` from the
            input cube label and use that DEM if present, otherwise
            fall back to ellipsoid. Matches ISIS cam2map's default.
          - ``"ellipsoid"`` or ``None``: use the constant mean radius
            from the cube's target body.
          - ``Path``: explicit path to an ISIS DEM cube.
    spice_source : {"isis", "naif", "auto"}
        Where ALE should source SPICE pointing/position data from when
        building the camera model:

          - ``"isis"`` (default): read the SPICE blobs **embedded in
            the cube**. This is the right choice for any pipeline that
            runs ``jigsaw update=true``, because jigsaw updates the
            cube's blobs but NOT the live NAIF kernels — reading the
            live kernels would silently throw away the bundle
            adjustment results.
          - ``"naif"``: force live NAIF kernel reads. Use only when
            comparing against pre-jigsaw geometry.
          - ``"auto"``: let ALE pick (currently prefers NAIF). Not
            recommended for jigsaw pipelines.
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

    timings: dict[str, float] | None = {} if profile else None
    t_total0 = time.perf_counter() if profile else None

    # Step 1: Load camera model + target body info. Default
    # spice_source="isis" reads the cube's embedded SPICE blobs, which is
    # the only correct choice after jigsaw update=true.
    console.print(
        f"[bold]Loading CSM camera model[/bold] from {input_cube.name} "
        f"(SPICE source: {spice_source})"
    )
    with _stage(timings, "load_camera"):
        model, body = load_camera(input_cube, spice_source=spice_source)
    n_lines, n_samples = get_image_size(model)
    console.print(f"  Input image: {n_samples} x {n_lines} (samples x lines)")
    console.print(
        f"  Target: {body.name} (NAIF {body.naif_id})  "
        f"radii eq={body.radius_equatorial_m:.1f} m  polar={body.radius_polar_m:.1f} m"
    )

    # Step 2: Resolve shape model. The DEM sampler's per-pixel radius
    # lookup falls back to the body's mean radius wherever the DEM
    # reports nodata. When the user disables the DEM with shape_model=
    # "ellipsoid", the same mean radius is the only surface radius used.
    mean_radius = body.radius_mean_m
    with _stage(timings, "dem_open"):
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
    with _stage(timings, "build_grid"):
        grid = _build_grid(
            model=model,
            body=body,
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
        with _stage(timings, "coord_transform"):
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
        with _stage(timings, "coord_transform"):
            coord_map = compute_transform_coarse(
                model,
                grid,
                mean_radius,
                step=coarse_step,
                input_n_lines=n_lines,
                input_n_samples=n_samples,
                dem_sampler=dem_sampler,
            )

    # Step 4b: Optional extra footprint-polygon mask. NOTE: this is NOT
    # what ISIS cam2map does — cam2map ignores the polygon entirely.
    # The flag is kept as an escape hatch for downstream callers that
    # want a polygon-masked output for their own reasons.
    if clip_to_footprint:
        console.print(
            "[bold]Applying extra footprint-polygon mask[/bold] "
            "(note: cam2map does not use the polygon)"
        )
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
    with _stage(timings, "read_input"):
        data, _label = read_isis_cube_raw(input_cube)

    # Step 6: Resample
    console.print(f"[bold]Resampling[/bold] ({interpolation.value})")
    with _stage(timings, "resample"):
        projected = resample(data, coord_map, interpolation=interpolation, fill_value=np.nan)

    # Step 7: Write output
    console.print(f"[bold]Writing[/bold] {output_path.name}")
    with _stage(timings, "write_output"):
        if output_format == "geotiff":
            result = write_geotiff(output_path, projected, grid, nodata=0.0)
            pvl_path = write_mapping_pvl(output_path, grid, body)
            console.print(f"  Mapping sidecar: {pvl_path.name}")
        else:
            msg = "ISIS cube output not yet implemented"
            raise NotImplementedError(msg)

    if timings is not None:
        total = time.perf_counter() - t_total0
        console.print()
        console.print("[bold cyan]Stage timings[/bold cyan]")
        for k, v in timings.items():
            pct = 100 * v / total if total > 0 else 0
            console.print(f"  {k:20s}  {v:7.2f} s  ({pct:5.1f}%)")
        accounted = sum(timings.values())
        other = total - accounted
        console.print(f"  {'(other)':20s}  {other:7.2f} s  ({100 * other / total:5.1f}%)")
        console.print(f"  {'total':20s}  {total:7.2f} s")

    console.print(f"[green bold]Done![/green bold] -> {result}")
    return result


def _rasterize_footprint(input_cube: Path, grid: OutputGrid) -> np.ndarray:
    """Rasterize the cube's footprint polygon onto the output grid.

    Reads the footprint polygon stored in the cube (by ``footprintinit``)
    and returns a boolean mask with True inside the polygon.

    Note: this mask is NOT used to reproduce ISIS ``cam2map`` behavior.
    Empirical testing (see ``docs/csm2map-design.md``) showed that
    ``cam2map`` ignores the polygon entirely. This helper is only invoked
    when the user explicitly passes ``--clip-to-footprint`` as an escape
    hatch for downstream polygon-masked outputs.

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
    body: TargetBody,
    input_cube: Path,
    map_file: Path | None,
    projection: str | None,
    resolution: float | None,
    lat_range: tuple[float, float] | None,
    lon_range: tuple[float, float] | None,
) -> OutputGrid:
    """Build output grid from MAP file or explicit parameters.

    When ``map_file`` is omitted and ``projection`` is also omitted, the
    default projection is an equirectangular CRS built from the target
    body's own radii (``body.radius_equatorial_m`` / ``.radius_polar_m``).
    This is what makes csm2map body-agnostic: Mars gets Mars radii, the
    Moon gets Moon radii, Europa gets Europa radii — all pulled from
    ALE's ISD via :class:`TargetBody`, never from a hardcoded default.
    """
    if map_file is not None:
        return grid_from_map_file(
            map_file,
            camera_lat_range=lat_range,
            camera_lon_range=lon_range,
            resolution_override=resolution,
        )

    if lat_range is None or lon_range is None:
        # Try to derive from the camera model image corners
        lat_range, lon_range = _derive_ground_range(model)

    if projection is None:
        # Default to equirectangular centered on the image. Using the
        # image's center latitude as lat_ts makes the local scale
        # isotropic at the image center (pixels are square there);
        # using the image's center longitude as lon_0 keeps the map
        # coordinates small and numerically well-conditioned. This
        # matches what ISIS cam2map does when auto-generating the
        # projection from the camera model.
        center_lat = (lat_range[0] + lat_range[1]) / 2.0
        center_lon = (lon_range[0] + lon_range[1]) / 2.0
        projection = (
            f"+proj=eqc +lat_ts={center_lat} +lon_0={center_lon} "
            f"+a={body.radius_equatorial_m} +b={body.radius_polar_m} "
            f"+units=m +no_defs +type=crs"
        )

    if resolution is None:
        # Auto-compute from the camera model's native ground sample
        # distance, matching ISIS cam2map's default behavior when no
        # PixelResolution is specified in the MAP file.
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


def _derive_ground_range(model) -> tuple[tuple[float, float], tuple[float, float]]:
    """Derive lat/lon ground range by probing CSM model at image corners + edges.

    Returns ``(lat_min, lat_max), (lon_min, lon_max)`` in degrees
    (planetocentric, positive east).

    Handles the **antimeridian crossing** case correctly: MRO is a polar
    orbiter, so CTX routinely images strips that cross the ±180° boundary.
    A naive ``min(lons)``/``max(lons)`` on the raw ``arctan2`` output
    would produce a 358° range instead of a 2° range. We avoid this by
    computing the *circular mean* of the longitudes, measuring all
    offsets from that mean in wrapped [-180, +180] space, and deriving
    the bounds from the offset extremes. This works for any image
    narrower than 180° — including antimeridian crossings, prime-meridian
    crossings, and polar images where the longitude spread is large but
    the physical strip is narrow.
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

    # Also sample along the edges at ~100 points per side
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

    lats_arr = np.array(lats)
    lons_arr = np.array(lons)

    # Latitude is straightforward — no wraparound.
    lat_span = lats_arr.max() - lats_arr.min()
    buf_lat = lat_span * 0.01
    lat_min = float(lats_arr.min() - buf_lat)
    lat_max = float(lats_arr.max() + buf_lat)

    # Longitude requires circular statistics to handle the antimeridian.
    # arctan2 returns values in [-180, +180]. A CTX strip crossing the
    # ±180° boundary will have probes at e.g. +179° and -179°, giving
    # a naive span of 358° instead of the correct 2°.
    #
    # Fix: compute the circular mean (the "center direction" of all the
    # probe longitudes), then measure every probe as an offset from that
    # center in wrapped [-180, +180] space. The offset min/max give the
    # true tight bounding box no matter where the data sits on the globe.
    lon_rad = np.radians(lons_arr)
    center_lon = np.degrees(np.arctan2(np.mean(np.sin(lon_rad)), np.mean(np.cos(lon_rad))))

    # Offset each probe from the circular mean, wrapped to [-180, +180].
    offsets = (lons_arr - center_lon + 180.0) % 360.0 - 180.0
    buf_lon = (offsets.max() - offsets.min()) * 0.01
    lon_min = float(center_lon + offsets.min() - buf_lon)
    lon_max = float(center_lon + offsets.max() + buf_lon)

    return (
        (lat_min, lat_max),
        (lon_min, lon_max),
    )
