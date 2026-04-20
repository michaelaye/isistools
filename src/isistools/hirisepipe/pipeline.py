"""Full HiRISE RED processing pipeline — ISIS-free.

Pipeline order (matching ISIS):

Per CCD:
    ingest → hical → histitch → cubenorm → map-project

Mosaic:
    equalizer → automos (rasterio.merge on projected GeoTIFFs)

Each CCD is projected with its own camera model so that the mosaic
assembly operates in map coordinates — no pixel offset calculations.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np

from isistools.hirisepipe.cubenorm import cubenorm
from isistools.hirisepipe.hical import hical, hical_from_edr
from isistools.hirisepipe.stitch import stitch_channels

try:
    from rich.console import Console

    _console = Console(stderr=True)

    def _log(msg: str) -> None:
        _console.print(msg)
except ImportError:

    def _log(msg: str) -> None:
        print(msg, file=sys.stderr, flush=True)


def calibrate_ccd(
    channel0_path: str | Path,
    channel1_path: str | Path,
    *,
    output_path: str | Path | None = None,
    units: str = "DN",
    matrices_dir: Path | None = None,
    balance: bool = True,
    normalize: bool = True,
) -> np.ndarray:
    """Calibrate a single HiRISE CCD: hical both channels, stitch, cubenorm.

    Accepts PDS EDR files (.IMG) or ISIS cubes (.cub).

    Parameters
    ----------
    channel0_path, channel1_path : path-like
        Paths to channel 0 and 1 files.
    output_path : path-like, optional
        If provided, write the calibrated CCD as a TIFF.
    units : str
        Calibration units ("DN", "DN/US", or "IOF").
    matrices_dir : Path, optional
        Calibration matrices directory.
    balance : bool
        Apply balance correction during stitching.
    normalize : bool
        Apply column normalization after stitching.

    Returns
    -------
    np.ndarray
        Calibrated, stitched, normalized CCD image.
    """
    ch0 = Path(channel0_path)
    ch1 = Path(channel1_path)

    cal_func = hical_from_edr if ch0.suffix.lower() == ".img" else hical
    cal_kwargs = {"units": units, "matrices_dir": matrices_dir}

    _log(f"Calibrating {ch0.name}...")
    t0 = time.perf_counter()
    cal0 = cal_func(ch0, **cal_kwargs)
    _log(f"  Channel 0 in {time.perf_counter() - t0:.1f}s")

    _log(f"Calibrating {ch1.name}...")
    t0 = time.perf_counter()
    cal1 = cal_func(ch1, **cal_kwargs)
    _log(f"  Channel 1 in {time.perf_counter() - t0:.1f}s")

    _log("Stitching channels...")
    t0 = time.perf_counter()
    stitched = stitch_channels(cal0, cal1, balance=balance)
    _log(f"  Stitched {stitched.shape[1]}x{stitched.shape[0]} in {time.perf_counter() - t0:.1f}s")

    if normalize:
        _log("Normalizing columns...")
        t0 = time.perf_counter()
        stitched = cubenorm(stitched)
        _log(f"  Normalized in {time.perf_counter() - t0:.1f}s")

    if output_path is not None:
        _write_tiff(stitched, Path(output_path))

    return stitched


def project_ccd(
    calibrated: np.ndarray,
    geometry_source: str | Path,
    output_path: str | Path,
    *,
    map_file: str | Path | None = None,
    projection: str | None = None,
    resolution: float | None = None,
    coarse_step: int = 32,
    interpolation: str = "bicubic",
) -> Path:
    """Map-project a single calibrated CCD image to GeoTIFF.

    Parameters
    ----------
    calibrated : np.ndarray
        Calibrated CCD image from ``calibrate_ccd()``.
    geometry_source : path-like
        Camera geometry source — spiceinit'd ISIS cube (.cub) or
        PDS EDR (.IMG) if ALE has a PDS3 driver for the instrument.
    output_path : path-like
        Output GeoTIFF path.
    map_file, projection, resolution, coarse_step, interpolation :
        See ``ctxpipe.pipeline.ctx_project`` for details.

    Returns
    -------
    Path to the output GeoTIFF.
    """
    # Reuse ctx_project's infrastructure — it handles camera loading,
    # grid building, transform, resample, and TIFF writing.
    from isistools.ctxpipe.pipeline import ctx_project

    return ctx_project(
        calibrated,
        geometry_source,
        output_path,
        map_file=map_file,
        projection=projection,
        resolution=resolution,
        coarse_step=coarse_step,
        interpolation=interpolation,
    )


def _find_edr_channels(
    obsid: str,
    ccdno: int,
    search_dirs: list[Path] | None = None,
) -> tuple[Path, Path]:
    """Find channel 0 and 1 EDR files for a given CCD."""
    if search_dirs is None:
        search_dirs = []
        try:
            from planetarypy.config import config

            root = Path(config["storage_root"])
            search_dirs.append(root / "mro" / "hirise" / "edr" / obsid)
            search_dirs.append(root / "mro" / "pds" / obsid)
        except Exception:
            pass
        search_dirs.append(Path.cwd())

    for d in search_dirs:
        ch0 = d / f"{obsid}_RED{ccdno}_0.IMG"
        ch1 = d / f"{obsid}_RED{ccdno}_1.IMG"
        if ch0.exists() and ch1.exists():
            return ch0, ch1

    searched = ", ".join(str(d) for d in search_dirs)
    raise FileNotFoundError(
        f"EDR files for {obsid}_RED{ccdno} channels 0+1 not found. "
        f"Searched: {searched}. "
        f"Download with: plp hiedr {obsid} --ccds {ccdno}"
    )


def _find_geometry_source(
    obsid: str,
    ccdno: int,
    search_dirs: list[Path] | None = None,
) -> Path:
    """Find a camera geometry source (spiceinit'd cube or EDR) for a CCD.

    Prefers channel 0 cube, falls back to EDR.
    """
    if search_dirs is None:
        search_dirs = []
        try:
            from planetarypy.config import config

            root = Path(config["storage_root"])
            search_dirs.append(root / "mro" / "hirise" / "edr" / obsid)
            search_dirs.append(root / "mro" / "pds" / obsid)
        except Exception:
            pass
        search_dirs.append(Path.cwd())
        search_dirs.append(Path("/tmp"))

    # Prefer spiceinit'd cube (has embedded SPICE tables)
    for d in search_dirs:
        cube = d / f"{obsid}_RED{ccdno}_0.cub"
        if cube.exists():
            return cube

    # Fall back to EDR (needs ALE PDS3 driver or web SpiceQL)
    for d in search_dirs:
        edr = d / f"{obsid}_RED{ccdno}_0.IMG"
        if edr.exists():
            return edr

    raise FileNotFoundError(
        f"No geometry source for {obsid}_RED{ccdno}. Need a spiceinit'd .cub or .IMG file."
    )


def _estimate_ccd_memory(edr_path: Path) -> float:
    """Estimate peak memory per CCD worker in bytes.

    Each worker holds: raw image (int16) + 2 calibrated channels (float32)
    + stitched (float32) + cubenorm copy (float32).
    """
    import pvl

    label = pvl.load(str(edr_path))
    lines = int(label["IMAGE"]["LINES"])
    samples = int(label["IMAGE"]["LINE_SAMPLES"])
    # Peak: ~4 float32 arrays of (lines × 2*samples) during stitch+norm
    return lines * samples * 2 * 4 * 4  # ~4 arrays, 2x width, 4 bytes each


def _smart_max_workers(
    sample_edr: Path,
    n_ccds: int,
    memory_fraction: float = 0.8,
) -> int:
    """Calculate how many parallel workers fit in available memory.

    Uses 80% of currently available RAM by default.

    Returns
    -------
    int
        Number of workers (at least 1, at most n_ccds).
    """
    try:
        import psutil

        available = psutil.virtual_memory().available
    except ImportError:
        import os

        # Fallback: assume 8 GB available if psutil not installed
        available = 8 * 1024**3
        _log("[yellow]psutil not installed — assuming 8 GB available[/yellow]")

    budget = available * memory_fraction
    per_ccd = _estimate_ccd_memory(sample_edr)
    workers = max(1, int(budget / per_ccd))
    workers = min(workers, n_ccds, os.cpu_count() or 4)

    _log(
        f"  Memory: {available / 1e9:.1f} GB available, "
        f"~{per_ccd / 1e9:.1f} GB/CCD → {workers} workers"
    )
    return workers


def _calibrate_ccd_to_file(args: tuple) -> Path:
    """Worker for parallel calibrate_all. Returns output path."""
    ccdno, ch0, ch1, out_path, units, matrices_dir = args
    _log(f"  [RED{ccdno}] Starting...")
    t0 = time.perf_counter()
    calibrate_ccd(
        ch0,
        ch1,
        output_path=out_path,
        units=units,
        matrices_dir=matrices_dir,
    )
    _log(f"  [RED{ccdno}] Done in {time.perf_counter() - t0:.1f}s → {out_path.name}")
    return out_path


def calibrate_all(
    obsid: str,
    ccds: list[int] | None = None,
    *,
    output_dir: str | Path | None = None,
    units: str = "DN",
    matrices_dir: Path | None = None,
    search_dirs: list[Path] | None = None,
    parallel: bool = True,
    max_workers: int | None = None,
) -> list[Path]:
    """Calibrate all CCDs and write per-CCD TIFFs.

    This is the recommended workflow for large images: calibrate each
    CCD independently to a TIFF, then project each one separately
    with csm2map.  Avoids holding multiple CCDs in memory at once.

    Parameters
    ----------
    obsid : str
        HiRISE observation ID.
    ccds : list of int, optional
        RED CCD numbers.  Default [4, 5].
    output_dir : path-like, optional
        Directory for output TIFFs.  Defaults to current directory.
    units : str
        Calibration units.
    matrices_dir : Path, optional
        Calibration matrices directory.
    search_dirs : list of Path, optional
        Directories to search for EDR files.
    parallel : bool
        If True (default), calibrate CCDs in parallel.
    max_workers : int or "smart", optional
        Maximum parallel workers.  If None or ``"smart"`` (default),
        auto-calculates based on available memory (uses 80% of free RAM).
        Set explicitly to control parallelism.

    Returns
    -------
    list of Path
        Paths to the calibrated per-CCD TIFFs.
    """
    if ccds is None:
        ccds = [4, 5]
    if output_dir is None:
        output_dir = Path.cwd()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve EDR paths and build worker args
    ccd_args = []
    for ccdno in sorted(ccds):
        ch0, ch1 = _find_edr_channels(obsid, ccdno, search_dirs)
        out_path = output_dir / f"{obsid}_RED{ccdno}.cal.norm.tif"
        ccd_args.append((ccdno, ch0, ch1, out_path, units, matrices_dir))

    if parallel and len(ccd_args) > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        if max_workers is None or max_workers == "smart":
            max_workers = _smart_max_workers(
                ccd_args[0][1],
                len(ccd_args),  # first EDR path, n_ccds
            )

        _log(f"Calibrating {obsid} {len(ccd_args)} CCDs → {output_dir} ({max_workers} workers)")
        t0 = time.perf_counter()

        try:
            from rich.progress import (
                BarColumn,
                Progress,
                SpinnerColumn,
                TextColumn,
                TimeElapsedColumn,
            )

            with Progress(
                SpinnerColumn(),
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                TextColumn("{task.completed}/{task.total} CCDs"),
                TimeElapsedColumn(),
                console=_console,
            ) as progress:
                task = progress.add_task(
                    f"Calibrating ({max_workers} workers)",
                    total=len(ccd_args),
                )
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    futures = {executor.submit(_calibrate_ccd_to_file, a): a[0] for a in ccd_args}
                    output_paths = []
                    for future in as_completed(futures):
                        output_paths.append(future.result())
                        progress.update(task, advance=1)
        except ImportError:
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                output_paths = list(executor.map(_calibrate_ccd_to_file, ccd_args))

        output_paths.sort()
        _log(f"[green]All CCDs calibrated in {time.perf_counter() - t0:.1f}s[/green]")
    else:
        _log(f"Calibrating {obsid} {len(ccd_args)} CCDs → {output_dir}")
        output_paths = []
        for args in ccd_args:
            output_paths.append(_calibrate_ccd_to_file(args))

    return output_paths


def create_red_mosaic(
    obsid: str,
    ccds: list[int] | None = None,
    *,
    units: str = "DN",
    matrices_dir: Path | None = None,
    search_dirs: list[Path] | None = None,
    output_path: str | Path | None = None,
    project: bool = False,
    resolution: float | None = None,
    projection: str | None = None,
    parallel: bool = True,
    max_workers: int | None = None,
) -> np.ndarray | Path:
    """Create a HiRISE RED mosaic — fully ISIS-free.

    Parameters
    ----------
    obsid : str
        HiRISE observation ID, e.g. "ESP_021491_0950".
    ccds : list of int, optional
        RED CCD numbers.  Default [4, 5] (central pair).
    units : str
        Calibration units.
    matrices_dir : Path, optional
        Calibration matrices directory.
    search_dirs : list of Path, optional
        Directories to search for EDR/cube files.
    output_path : path-like, optional
        Output file path.  For unprojected: writes a TIFF.
        For projected: writes a GeoTIFF mosaic.
    project : bool
        If True, map-project each CCD before mosaicking.
        Requires a spiceinit'd cube per CCD for camera geometry.
    resolution : float, optional
        Pixel resolution in m/px for projection (auto if omitted).
    projection : str, optional
        PROJ string for projection (auto-selected if omitted).
    parallel : bool
        If True (default), calibrate CCDs in parallel using multiple
        CPU cores.  Set to False for sequential processing.
    max_workers : int, optional
        Maximum number of parallel workers.  Defaults to the number
        of CCDs or the CPU count, whichever is smaller.

    Returns
    -------
    np.ndarray (unprojected) or Path (projected)
    """
    if ccds is None:
        ccds = [4, 5]

    _log(f"HiRISE RED mosaic: {obsid}, CCDs {ccds}, projected={project}")

    if not project:
        return _create_raw_mosaic(
            obsid,
            ccds,
            units=units,
            matrices_dir=matrices_dir,
            search_dirs=search_dirs,
            output_path=output_path,
            parallel=parallel,
            max_workers=max_workers,
        )
    else:
        return _create_projected_mosaic(
            obsid,
            ccds,
            units=units,
            matrices_dir=matrices_dir,
            search_dirs=search_dirs,
            output_path=output_path,
            resolution=resolution,
            projection=projection,
        )


def _calibrate_one_ccd(args: tuple) -> tuple[int, np.ndarray]:
    """Worker function for parallel CCD calibration.

    Takes a tuple to be compatible with ProcessPoolExecutor.map().
    """
    ccdno, ch0_path, ch1_path, units, matrices_dir = args
    _log(f"  [RED{ccdno}] Starting calibration...")
    t0 = time.perf_counter()
    result = calibrate_ccd(ch0_path, ch1_path, units=units, matrices_dir=matrices_dir)
    _log(
        f"  [RED{ccdno}] Done in {time.perf_counter() - t0:.1f}s "
        f"({result.shape[1]}x{result.shape[0]})"
    )
    return (ccdno, result)


def _create_raw_mosaic(
    obsid: str,
    ccds: list[int],
    *,
    units: str,
    matrices_dir: Path | None,
    search_dirs: list[Path] | None,
    output_path: str | Path | None,
    parallel: bool = True,
    max_workers: int | None = None,
) -> np.ndarray:
    """Calibrate and mosaic in raw camera geometry (no projection)."""
    # Resolve EDR paths first (sequential, fast)
    ccd_args = []
    for ccdno in sorted(ccds):
        ch0, ch1 = _find_edr_channels(obsid, ccdno, search_dirs)
        ccd_args.append((ccdno, ch0, ch1, units, matrices_dir))

    # Calibrate CCDs
    if parallel and len(ccd_args) > 1:
        from concurrent.futures import ProcessPoolExecutor, as_completed

        if max_workers is None or max_workers == "smart":
            max_workers = _smart_max_workers(
                ccd_args[0][1],
                len(ccd_args),  # first EDR path, n_ccds
            )

        t0 = time.perf_counter()
        try:
            from rich.progress import (
                BarColumn,
                Progress,
                SpinnerColumn,
                TextColumn,
                TimeElapsedColumn,
            )

            with Progress(
                SpinnerColumn(),
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                TextColumn("{task.completed}/{task.total} CCDs"),
                TimeElapsedColumn(),
                console=_console,
            ) as progress:
                task = progress.add_task(
                    f"Calibrating ({max_workers} workers)",
                    total=len(ccd_args),
                )
                with ProcessPoolExecutor(max_workers=max_workers) as executor:
                    futures = {executor.submit(_calibrate_one_ccd, a): a[0] for a in ccd_args}
                    results = []
                    for future in as_completed(futures):
                        ccdno, img = future.result()
                        results.append((ccdno, img))
                        progress.update(task, advance=1)
        except ImportError:
            _log(f"\nCalibrating {len(ccd_args)} CCDs in parallel ({max_workers} workers)...")
            with ProcessPoolExecutor(max_workers=max_workers) as executor:
                results = list(executor.map(_calibrate_one_ccd, ccd_args))

        results.sort(key=lambda x: x[0])
        ccd_images = [img for _, img in results]
        _log(f"[green]All CCDs calibrated in {time.perf_counter() - t0:.1f}s[/green]")
    else:
        ccd_images = []
        for args in ccd_args:
            _, img = _calibrate_one_ccd(args)
            ccd_images.append(img)

    if len(ccd_images) == 1:
        mosaic = ccd_images[0]
    else:
        _log("\nAssembling raw mosaic (48px overlap)...")
        mosaic = _assemble_raw(ccd_images)

    if output_path is not None:
        _write_tiff(mosaic, Path(output_path))

    return mosaic


def _create_projected_mosaic(
    obsid: str,
    ccds: list[int],
    *,
    units: str,
    matrices_dir: Path | None,
    search_dirs: list[Path] | None,
    output_path: str | Path | None,
    resolution: float | None,
    projection: str | None,
) -> Path:
    """Calibrate, project each CCD, then mosaic in map coordinates."""
    import tempfile

    projected_paths = []
    tmpdir = Path(tempfile.mkdtemp(prefix="hirisepipe_"))

    for ccdno in sorted(ccds):
        _log(f"\n--- CCD RED{ccdno} ---")

        # Calibrate
        ch0, ch1 = _find_edr_channels(obsid, ccdno, search_dirs)
        cal = calibrate_ccd(ch0, ch1, units=units, matrices_dir=matrices_dir)

        # Find camera geometry for this specific CCD
        geom = _find_geometry_source(obsid, ccdno, search_dirs)
        _log(f"Geometry: {geom.name}")

        # Project this CCD
        ccd_tif = tmpdir / f"{obsid}_RED{ccdno}_projected.tif"
        project_ccd(
            cal,
            geom,
            ccd_tif,
            resolution=resolution,
            projection=projection,
        )
        projected_paths.append(ccd_tif)

    # Mosaic the projected GeoTIFFs
    if output_path is None:
        output_path = Path.cwd() / f"{obsid}_RED{''.join(str(c) for c in ccds)}.tif"
    output_path = Path(output_path)

    if len(projected_paths) == 1:
        import shutil

        shutil.move(str(projected_paths[0]), str(output_path))
    else:
        _log("\nMosaicking projected CCDs...")
        _merge_geotiffs(projected_paths, output_path)

    # Clean up temp files
    for p in projected_paths:
        p.unlink(missing_ok=True)
    tmpdir.rmdir()

    _log(f"\nMosaic: {output_path}")
    return output_path


def _merge_geotiffs(input_paths: list[Path], output_path: Path) -> None:
    """Merge multiple GeoTIFFs into one using rasterio."""
    import rasterio
    from rasterio.merge import merge

    datasets = [rasterio.open(str(p)) for p in input_paths]
    try:
        mosaic, transform = merge(datasets)
        profile = datasets[0].profile.copy()
        profile.update(
            height=mosaic.shape[1],
            width=mosaic.shape[2],
            transform=transform,
            compress="zstd",
        )
        with rasterio.open(str(output_path), "w", **profile) as dst:
            dst.write(mosaic)
    finally:
        for ds in datasets:
            ds.close()


def _assemble_raw(
    ccd_images: list[np.ndarray],
    overlap_px: int = 48,
) -> np.ndarray:
    """Assemble CCD images in raw geometry with pixel overlap."""
    if len(ccd_images) == 1:
        return ccd_images[0]

    max_lines = max(img.shape[0] for img in ccd_images)
    widths = [img.shape[1] for img in ccd_images]
    total_width = sum(widths) - overlap_px * (len(ccd_images) - 1)

    mosaic = np.full((max_lines, total_width), np.nan, dtype=np.float32)
    col = 0
    for img in ccd_images:
        h, w = img.shape
        mosaic[:h, col : col + w] = img
        col += w - overlap_px

    return mosaic


def _write_tiff(image: np.ndarray, output_path: Path) -> None:
    """Write a float32 TIFF."""
    import rasterio

    _log(f"Writing {output_path.name}...")
    height, width = image.shape
    with rasterio.open(
        str(output_path),
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        compress="zstd",
        tiled=True,
        blockxsize=256,
        blockysize=256,
    ) as dst:
        dst.write(np.nan_to_num(image, nan=0.0), 1)
        dst.nodata = 0.0
