"""Full CTX processing pipeline — ingest, calibrate, even/odd correct, project.

Orchestrates the complete CTX processing chain as composable functions.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from isistools.ctxpipe.calibrate import calibrate
from isistools.ctxpipe.evenodd import correct_evenodd
from isistools.ctxpipe.ingest import CTXMetadata, ingest_ctx_edr


def ctx_calibrate(
    edr_path: str | Path,
    *,
    iof: bool = False,
    sun_distance_km: float | None = None,
    calibration_dir: Path | None = None,
    evenodd: bool = True,
    fill_gap: bool = True,
) -> tuple[np.ndarray, CTXMetadata]:
    """Calibrate a CTX EDR: ingest → dark/flat → optional even/odd correction.

    This is the main calibration entry point.  Returns a Level 1 calibrated
    image (unprojected, in camera geometry).

    Parameters
    ----------
    edr_path : path-like
        Path to the PDS3 CTX EDR (.IMG file).
    iof : bool
        Convert to I/F units.  Requires ``sun_distance_km`` or spiceypy.
    sun_distance_km : float, optional
        Sun-to-target distance in km.  Auto-computed from a co-located
        spiceinit'd cube if spiceypy is installed.
    calibration_dir : Path, optional
        Directory containing CTX calibration files.
    evenodd : bool
        Apply even/odd column correction (default True).  Only has an
        effect for spatial summing = 1.  Set to False for debugging or
        if downstream processing handles it.
    fill_gap : bool
        Treat raw DN 0 as NULL (NaN).

    Returns
    -------
    image : np.ndarray
        Calibrated float32 image.  NaN where data is missing.
    metadata : CTXMetadata
        Extracted metadata.
    """
    import sys
    import time

    def _log(msg: str) -> None:
        print(msg, file=sys.stderr, flush=True)

    # Step 1: Ingest
    _log(f"Ingesting {Path(edr_path).name}...")
    t0 = time.perf_counter()
    raw, metadata = ingest_ctx_edr(edr_path, fill_gap=fill_gap)
    _log(f"  {raw.shape[1]}x{raw.shape[0]} pixels in {time.perf_counter() - t0:.1f}s")

    # Auto-compute Sun distance for I/F
    if iof and sun_distance_km is None:
        _log("Computing Sun distance...")
        sun_distance_km = _auto_sun_distance(edr_path)

    # Step 2: Calibrate
    _log("Calibrating...")
    t0 = time.perf_counter()
    cal = calibrate(
        raw,
        metadata,
        iof=iof,
        sun_distance_km=sun_distance_km,
        calibration_dir=calibration_dir,
    )
    _log(f"  Calibrated in {time.perf_counter() - t0:.1f}s")

    # Step 3: Even/odd correction
    if evenodd:
        _log("Applying even/odd correction...")
        t0 = time.perf_counter()
        cal = correct_evenodd(cal, spatial_summing=metadata.spatial_summing)
        _log(f"  Even/odd in {time.perf_counter() - t0:.1f}s")

    return cal, metadata


def ctx_project(
    calibrated: np.ndarray,
    geometry_source: str | Path,
    output_path: str | Path,
    *,
    map_file: str | Path | None = None,
    projection: str | None = None,
    resolution: float | None = None,
    coarse_step: int = 32,
    interpolation: str = "bicubic",
    pvl_sidecar: bool = False,
) -> Path:
    """Map-project a calibrated CTX image to a GeoTIFF.

    The camera model is loaded from either a PDS EDR (via ALE + NAIF
    kernels, no ISIS needed) or a spiceinit'd ISIS cube (fallback).

    Default projection: Sinusoidal for |lat| < 70, Polar Stereographic
    for |lat| >= 70.

    Parameters
    ----------
    calibrated : np.ndarray
        Calibrated image from ``ctx_calibrate()``.
    geometry_source : path-like
        Source for camera geometry.  Can be:

        - A PDS3 EDR (.IMG) — ALE builds the camera model directly
          from the label + NAIF kernels.  No spiceinit needed.
        - A spiceinit'd ISIS cube (.cub) — ALE reads the embedded
          SPICE tables.  Pixel data is ignored.
        - A pre-computed ISD (.json) — loaded directly, no ALE call.
    output_path : path-like
        Output GeoTIFF path.
    map_file : path-like, optional
        ISIS MAP file — overrides automatic projection selection.
    projection : str, optional
        PROJ string — overrides automatic projection selection.
    resolution : float, optional
        Pixel resolution in m/px.  Auto-derived from camera GSD if omitted.
    coarse_step : int
        Coarse grid step for coordinate transform interpolation.
    interpolation : str
        Resampling: "nearest", "bilinear", or "bicubic".
    pvl_sidecar : bool
        Write ISIS-compatible PVL sidecar (default False).

    Returns
    -------
    Path
        Path to the output GeoTIFF.
    """
    import sys
    import time

    from isistools.csm2map.camera import (
        compute_ground_sample_distance,
        get_image_size,
    )
    from isistools.csm2map.grid import grid_from_map_file, grid_from_params
    from isistools.csm2map.pipeline import _derive_ground_range
    from isistools.csm2map.resample import Interpolation, resample
    from isistools.csm2map.transform import compute_transform_coarse
    from isistools.csm2map.writers import write_geotiff, write_mapping_pvl

    def _log(msg: str) -> None:
        print(msg, file=sys.stderr, flush=True)

    output_path = Path(output_path)

    # Load camera model — try direct label path first, fall back to cube
    _log("Loading camera model...")
    t0 = time.perf_counter()
    model, body = _load_camera_auto(geometry_source)
    n_lines, n_samples = get_image_size(model)
    _log(
        f"  Camera loaded in {time.perf_counter() - t0:.1f}s ({body.name}, {n_samples}x{n_lines})"
    )

    cal_lines, cal_samps = calibrated.shape
    if cal_lines != n_lines:
        raise ValueError(
            f"Calibrated image lines ({cal_lines}) does not match camera model lines ({n_lines})"
        )
    # EDR-path models report full line width (including prefix/suffix),
    # while the calibrated image has science pixels only.  The CSM model
    # coordinates are aligned to the science pixels (sample 0 = first
    # science pixel), so we use the calibrated image dimensions for the
    # transform, not the model's reported dimensions.
    n_samples = cal_samps

    mean_radius = body.radius_mean_m

    # Build output grid
    _log("Building output grid...")
    if map_file is not None:
        grid = grid_from_map_file(map_file)
    else:
        lat_range, lon_range = _derive_ground_range(model)

        if resolution is None:
            resolution = compute_ground_sample_distance(model, body)

        if projection is None:
            projection = _auto_projection(lat_range, lon_range, body)

        grid = grid_from_params(
            crs=projection,
            resolution=resolution,
            lat_min=lat_range[0],
            lat_max=lat_range[1],
            lon_min=lon_range[0],
            lon_max=lon_range[1],
        )
    _log(f"  {grid.width}x{grid.height} pixels, {grid.resolution:.2f} m/px")

    # Coordinate transform + resample
    _log("Computing coordinate transform...")
    t0 = time.perf_counter()
    coord_map = compute_transform_coarse(
        model,
        grid,
        mean_radius,
        step=coarse_step,
        input_n_lines=n_lines,
        input_n_samples=n_samples,
    )
    _log(f"  Transform in {time.perf_counter() - t0:.1f}s")

    _log("Resampling...")
    t0 = time.perf_counter()
    interp_enum = Interpolation(interpolation)
    projected = resample(calibrated, coord_map, interpolation=interp_enum, fill_value=np.nan)
    _log(f"  Resample in {time.perf_counter() - t0:.1f}s")

    # Write output
    _log("Writing GeoTIFF...")
    result = write_geotiff(output_path, projected, grid, nodata=0.0)
    if pvl_sidecar:
        write_mapping_pvl(output_path, grid, body)
    return result


def ctx_edr_to_map(
    edr_path: str | Path,
    output_path: str | Path,
    *,
    geometry_source: str | Path | None = None,
    iof: bool = False,
    sun_distance_km: float | None = None,
    calibration_dir: Path | None = None,
    evenodd: bool = True,
    map_file: str | Path | None = None,
    projection: str | None = None,
    resolution: float | None = None,
    interpolation: str = "bicubic",
    pvl_sidecar: bool = False,
) -> Path:
    """Full pipeline: PDS EDR → calibrated → map-projected GeoTIFF.

    Combines ``ctx_calibrate`` and ``ctx_project`` in one call.
    By default uses the EDR itself for camera geometry (via ALE +
    NAIF kernels).  No spiceinit or ISIS cube needed.

    Parameters
    ----------
    edr_path : path-like
        PDS3 CTX EDR (.IMG file).
    output_path : path-like
        Output GeoTIFF path.
    geometry_source : path-like, optional
        Source for camera geometry.  Defaults to the EDR itself (ALE
        builds the camera model from NAIF kernels).  Can also be a
        spiceinit'd ISIS cube or a pre-computed ISD JSON.
    iof : bool
        Convert to I/F units.
    sun_distance_km : float, optional
        Sun-target distance for I/F.
    calibration_dir : Path, optional
        CTX calibration data directory.
    evenodd : bool
        Apply even/odd column correction (default True).
    map_file : path-like, optional
        ISIS MAP file for projection.
    projection : str, optional
        PROJ string for projection.
    resolution : float, optional
        Pixel resolution in m/px.
    interpolation : str
        Resampling method.
    pvl_sidecar : bool
        Write ISIS-compatible PVL sidecar (default False).

    Returns
    -------
    Path
        Path to the output GeoTIFF.
    """
    if geometry_source is None:
        geometry_source = edr_path

    cal, _meta = ctx_calibrate(
        edr_path,
        iof=iof,
        sun_distance_km=sun_distance_km,
        calibration_dir=calibration_dir,
        evenodd=evenodd,
    )
    return ctx_project(
        cal,
        geometry_source,
        output_path,
        map_file=map_file,
        projection=projection,
        resolution=resolution,
        interpolation=interpolation,
        pvl_sidecar=pvl_sidecar,
    )


def _load_camera_auto(geometry_source: str | Path):
    """Load a CSM camera model from an EDR, cube, or cached ISD.

    Tries in order:
    1. load_camera_from_label (NaifSpice — works with raw EDRs)
    2. load_camera (IsisSpice — needs spiceinit'd cube)
    """
    from pathlib import Path

    geometry_source = Path(geometry_source)
    suffix = geometry_source.suffix.lower()

    # If it's a .json ISD, load directly
    if suffix == ".json":
        from isistools.csm2map.camera import load_camera_from_label

        return load_camera_from_label(geometry_source, refresh_isd=False)

    # Try NaifSpice first (works with EDRs, no spiceinit needed)
    try:
        from isistools.csm2map.camera import load_camera_from_label

        return load_camera_from_label(geometry_source)
    except Exception:
        pass

    # Fall back to IsisSpice (needs spiceinit'd cube)
    try:
        from isistools.csm2map.camera import load_camera

        return load_camera(geometry_source)
    except Exception as e:
        raise RuntimeError(
            f"Could not load camera model from {geometry_source}. "
            f"Provide a PDS EDR (with NAIF kernels available) or a "
            f"spiceinit'd ISIS cube. Error: {e}"
        ) from e


def _auto_projection(
    lat_range: tuple[float, float],
    lon_range: tuple[float, float],
    body,
) -> str:
    """Select projection based on latitude range, centered on the image.

    - |lat| < 70: Sinusoidal (area-preserving, good for mid-latitudes)
    - |lat| >= 70: Polar Stereographic (avoids singularity at poles)

    The projection's central meridian (``lon_0``) is set to the mean
    image longitude.  This keeps the image's axis-aligned bounding
    box in projected coordinates tight against the actual ground
    footprint.  Leaving ``lon_0 = 0`` (as earlier versions did) made
    any image far from the prime meridian explode its bounding box
    — e.g. a CTX strip at lon 102° with a 5° latitude span would
    occupy ~470 km in projected X even though the ground footprint
    was only ~30 km wide, which in turn made
    ``compute_transform_coarse`` allocate the transform arrays on a
    4 billion-pixel grid.
    """
    center_lat = (lat_range[0] + lat_range[1]) / 2.0
    center_lon = (lon_range[0] + lon_range[1]) / 2.0

    a = body.radius_equatorial_m
    b = body.radius_polar_m

    if abs(center_lat) >= 70:
        # Polar Stereographic
        lat_0 = 90.0 if center_lat > 0 else -90.0
        return (
            f"+proj=stere +lat_0={lat_0} +lon_0={center_lon} "
            f"+a={a} +b={b} +units=m +no_defs +type=crs"
        )
    else:
        # Sinusoidal (ISIS default)
        return (
            f"+proj=sinu +lon_0={center_lon} +a={a} +b={b} +units=m +no_defs +type=crs"
        )


def _auto_sun_distance(edr_path: str | Path) -> float:
    """Try to auto-compute Sun distance from a co-located spiceinit'd cube."""
    try:
        from isistools.spice_utils import sun_distance_from_cube
    except ImportError:
        raise ValueError(
            "iof=True requires sun_distance_km or spiceypy. "
            "Install spiceypy: conda install -c conda-forge spiceypy"
        )

    edr = Path(edr_path)
    cube_candidates = [
        edr.with_suffix(".cub"),
        edr.parent / (edr.stem + ".raw.cub"),
    ]
    for candidate in cube_candidates:
        if candidate.exists():
            return sun_distance_from_cube(candidate)

    raise ValueError(
        "iof=True requires sun_distance_km. Either provide it "
        "explicitly or ensure a spiceinit'd .cub file exists "
        "alongside the EDR."
    )
