"""Standalone CLI for ctxpipe.

Entry point: ``ctxpipe`` at the shell (registered in pyproject.toml).

Two commands:

- ``ctxpipe calibrate INPUT.IMG OUTPUT.tif``  — Level 1 calibrated TIFF
- ``ctxpipe project INPUT.IMG CUBE.cub OUTPUT.tif`` — Level 2 map-projected GeoTIFF
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

app = typer.Typer(
    name="ctxpipe",
    help="Python CTX calibration pipeline (replaces ISIS CTX processing chain).",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _write_preview_png(image, output_path: Path, plow: float = 1, phigh: float = 99):
    """Write a percentile-stretched 8-bit PNG for quick viewing."""
    import numpy as np
    from PIL import Image

    valid = image[np.isfinite(image)]
    if len(valid) == 0:
        return

    vmin = np.percentile(valid, plow)
    vmax = np.percentile(valid, phigh)
    if vmax == vmin:
        vmax = vmin + 1

    stretched = np.clip((image - vmin) / (vmax - vmin) * 255, 0, 255)
    stretched = np.nan_to_num(stretched, nan=0).astype(np.uint8)

    png_path = output_path.with_suffix(".png")
    Image.fromarray(stretched).save(png_path)
    return png_path


@app.command()
def calibrate(
    from_edr: Annotated[
        Path,
        typer.Argument(help="Input PDS3 CTX EDR file (.IMG)."),
    ],
    to: Annotated[
        Path,
        typer.Argument(help="Output calibrated TIFF file."),
    ],
    iof: Annotated[
        bool,
        typer.Option("--iof/--no-iof", help="Convert to I/F units."),
    ] = False,
    no_evenodd: Annotated[
        bool,
        typer.Option("--no-evenodd", help="Skip even/odd column correction."),
    ] = False,
    sun_distance: Annotated[
        Optional[float],
        typer.Option("--sun-distance", help="Sun-target distance in km (for --iof)."),
    ] = None,
    cal_dir: Annotated[
        Optional[Path],
        typer.Option("--cal-dir", help="CTX calibration data directory."),
    ] = None,
    png: Annotated[
        bool,
        typer.Option("--png", help="Write a stretched 8-bit PNG preview alongside."),
    ] = False,
    no_compress: Annotated[
        bool,
        typer.Option("--no-compress", help="Skip TIFF compression (faster write)."),
    ] = False,
) -> None:
    """Produce a Level 1 calibrated TIFF from a CTX EDR.

    Reads a PDS3 CTX EDR, applies SQROOT decompression, dark current
    subtraction, flat-field correction, and (by default) even/odd column
    correction.  Output is an unprojected TIFF in camera geometry with
    calibrated DN values (or I/F with --iof).

    Even/odd correction is applied automatically for summing=1 images
    (where the CTX detector exhibits alternating-column striping).
    Use --no-evenodd to skip it for debugging or custom downstream
    processing.

    The output TIFF is float32 (not viewable in macOS Preview). Use
    --png to write a percentile-stretched 8-bit PNG preview alongside.
    Use --no-compress for faster write at the cost of larger files.
    """
    import numpy as np
    import rasterio

    from isistools.ctxpipe.pipeline import ctx_calibrate

    image, meta = ctx_calibrate(
        from_edr,
        iof=iof,
        sun_distance_km=sun_distance,
        calibration_dir=cal_dir,
        evenodd=not no_evenodd,
    )

    to = Path(to)
    height, width = image.shape
    write_opts = {
        "driver": "GTiff",
        "height": height,
        "width": width,
        "count": 1,
        "dtype": "float32",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }
    if not no_compress:
        write_opts["compress"] = "zstd"

    with rasterio.open(str(to), "w", **write_opts) as dst:
        dst.write(np.nan_to_num(image, nan=0.0), 1)
        dst.update_tags(
            PRODUCT_ID=meta.product_id,
            TARGET=meta.target_name,
            EXPOSURE_MS=str(meta.line_exposure_duration),
            SPATIAL_SUMMING=str(meta.spatial_summing),
        )
        dst.nodata = 0.0

    valid = image[np.isfinite(image)]
    typer.echo(
        f"{meta.product_id}: {width}x{height}, "
        f"mean={valid.mean():.3f}, "
        f"evenodd={'on' if not no_evenodd else 'off'}"
    )

    if png:
        png_path = _write_preview_png(image, to)
        typer.echo(f"Preview: {png_path}")


@app.command()
def project(
    from_edr: Annotated[
        Path,
        typer.Argument(help="Input PDS3 CTX EDR file (.IMG)."),
    ],
    to: Annotated[
        Path,
        typer.Argument(help="Output map-projected GeoTIFF."),
    ],
    geometry: Annotated[
        Optional[Path],
        typer.Option(
            "--geometry",
            "-g",
            help="Camera geometry source: EDR, spiceinit'd cube, or ISD JSON. "
            "Defaults to the input EDR (uses ALE + NAIF kernels).",
        ),
    ] = None,
    iof: Annotated[
        bool,
        typer.Option("--iof/--no-iof", help="Convert to I/F units."),
    ] = False,
    no_evenodd: Annotated[
        bool,
        typer.Option("--no-evenodd", help="Skip even/odd column correction."),
    ] = False,
    resolution: Annotated[
        Optional[float],
        typer.Option("-r", "--resolution", help="Pixel resolution in m/px (auto if omitted)."),
    ] = None,
    map_file: Annotated[
        Optional[Path],
        typer.Option("--map", help="ISIS MAP file for projection definition."),
    ] = None,
    projection: Annotated[
        Optional[str],
        typer.Option("--projection", help="PROJ string (overrides auto-selection)."),
    ] = None,
    interpolation: Annotated[
        str,
        typer.Option("--interp", help="Resampling: nearest, bilinear, bicubic."),
    ] = "bicubic",
    pvl_sidecar: Annotated[
        bool,
        typer.Option("--pvl-sidecar", help="Write ISIS-compatible PVL sidecar."),
    ] = False,
    png: Annotated[
        bool,
        typer.Option("--png", help="Write a stretched 8-bit PNG preview alongside."),
    ] = False,
    tile_size: Annotated[
        str,
        typer.Option(
            "--tile-size",
            help="Output tiling: 'auto' (default; size by available RAM), "
            "'none' (force batch path), or a positive integer for a fixed "
            "tile edge in pixels. Auto falls through to batch when one "
            "tile covers the whole output.",
        ),
    ] = "auto",
) -> None:
    """Produce a Level 2 map-projected GeoTIFF from a CTX EDR.

    Calibrates the EDR, then map-projects using a CSM camera model.
    By default, ALE builds the camera model directly from the EDR
    label and NAIF kernels — no spiceinit or ISIS cube needed.

    Use -g/--geometry to specify an alternative geometry source
    (spiceinit'd cube for jigsaw-adjusted pointing, or a cached ISD).

    Default projection is auto-selected based on latitude:
    Sinusoidal for |lat| < 70, Polar Stereographic for |lat| >= 70.
    Override with --map (ISIS MAP file) or --projection (PROJ string).

    Resolution is auto-derived from the camera's ground sample distance
    unless explicitly set with -r/--resolution.

    The output GeoTIFF is float32. Use --png for a quick-look preview.
    """
    from isistools.ctxpipe.pipeline import ctx_calibrate, ctx_project

    cal, _meta = ctx_calibrate(from_edr, iof=iof, evenodd=not no_evenodd)

    geometry_source = geometry if geometry is not None else from_edr

    to = Path(to)
    result = ctx_project(
        cal,
        geometry_source,
        to,
        map_file=map_file,
        projection=projection,
        resolution=resolution,
        interpolation=interpolation,
        pvl_sidecar=pvl_sidecar,
        tile_size=tile_size,
    )
    typer.echo(f"Projected -> {result}")

    if png:
        import rasterio

        with rasterio.open(str(result)) as src:
            data = src.read(1)
        png_path = _write_preview_png(data, to)
        typer.echo(f"Preview: {png_path}")
