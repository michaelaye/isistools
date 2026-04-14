"""Standalone CLI for csm2map.

Entry point: ``csm2map`` at the shell (registered in pyproject.toml).
Two commands:

- ``csm2map input.cub output.tif [OPTIONS]``   — map-project
- ``csm2map compare isis.cub csm.tif``         — validate against ISIS
"""

from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer

app = typer.Typer(
    name="csm2map",
    help="CSM-based map projection for planetary images (ISIS cam2map replacement).",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
)


@app.command()
def project(
    from_cube: Annotated[
        Path,
        typer.Argument(help="Input spiceinit'd ISIS cube.", exists=True),
    ],
    to: Annotated[
        Path,
        typer.Argument(help="Output file path (.tif for GeoTIFF)."),
    ],
    map: Annotated[
        Optional[Path],
        typer.Option(
            "--map", help="ISIS MAP PVL file defining projection/resolution/range.", exists=True
        ),
    ] = None,
    projection: Annotated[
        Optional[str],
        typer.Option("--projection", help="PROJ string (default: equirectangular)."),
    ] = None,
    resolution: Annotated[
        Optional[float],
        typer.Option(
            "--resolution",
            "-r",
            help="Pixel resolution in meters/pixel. If omitted, auto-computed "
            "from the camera model's ground sample distance at the image center "
            "(matching ISIS cam2map's default behavior).",
        ),
    ] = None,
    minlat: Annotated[
        Optional[float],
        typer.Option("--minlat", help="Minimum latitude (degrees)."),
    ] = None,
    maxlat: Annotated[
        Optional[float],
        typer.Option("--maxlat", help="Maximum latitude (degrees)."),
    ] = None,
    minlon: Annotated[
        Optional[float],
        typer.Option("--minlon", help="Minimum longitude (degrees)."),
    ] = None,
    maxlon: Annotated[
        Optional[float],
        typer.Option("--maxlon", help="Maximum longitude (degrees)."),
    ] = None,
    step: Annotated[
        int,
        typer.Option(
            "--step",
            "-s",
            help="Coarse grid step in output pixels. Smaller=more accurate but "
            "slower (more CSM calls). Default 32 works well for CTX at 6 m/pix "
            "and any similarly smooth line-scan + DEM combination; drop to 16 "
            "or 8 for HiRISE / rugged terrain / high-res DEM.",
        ),
    ] = 32,
    dense: Annotated[
        bool,
        typer.Option("--dense", help="Evaluate CSM at every pixel (slow, for validation)."),
    ] = False,
    validate: Annotated[
        bool,
        typer.Option("--validate", help="Spot-check coarse transform accuracy."),
    ] = False,
    clip_to_footprint: Annotated[
        bool,
        typer.Option(
            "--clip-to-footprint",
            help="Apply an extra mask from the footprint polygon stored by "
            "footprintinit. This does NOT reproduce ISIS cam2map behavior — "
            "cam2map ignores the polygon entirely. Use only if you want a "
            "polygon-clipped output for your own downstream reasons.",
        ),
    ] = False,
    shape_model: Annotated[
        str,
        typer.Option(
            "--shape-model",
            help="Shape model: 'auto' (read from cube label, default), "
            "'ellipsoid' (constant mean radius), or path to a DEM cube.",
        ),
    ] = "auto",
    spice_source: Annotated[
        str,
        typer.Option(
            "--spice-source",
            help="SPICE pointing source: 'isis' (cube's embedded blobs, "
            "the only correct choice after jigsaw update=true), "
            "'naif' (live kernels), or 'auto' (ALE picks, prefers NAIF). "
            "Default 'isis'.",
        ),
    ] = "isis",
    profile: Annotated[
        bool,
        typer.Option("--profile", help="Print per-stage wall times at the end."),
    ] = False,
    interp: Annotated[
        str,
        typer.Option("--interp", "-i", help="Interpolation: nearest, bilinear, bicubic."),
    ] = "bicubic",
):
    """Map-project an ISIS cube using a CSM camera model.

    Drop-in replacement for ISIS cam2map. By default, reads the shape model
    from the cube label, auto-computes resolution from the camera's ground
    sample distance, derives bounds from the camera footprint, and centers
    the projection on the image. All five parameters are derived from the
    cube itself — the only required arguments are input and output paths.

    \b
    Examples:
      csm2map input.cub output.tif                    # auto everything
      csm2map input.cub output.tif --map equi.map     # use an ISIS MAP file
      csm2map input.cub output.tif -r 6.0             # explicit 6 m/px
      csm2map input.cub output.tif --dense --validate
      csm2map input.cub output.tif --shape-model ellipsoid
    """
    try:
        from isistools.csm2map.pipeline import csm2map as _csm2map
        from isistools.csm2map.resample import Interpolation
    except ImportError as e:
        typer.echo(f"csm2map requires extra dependencies: pip install isistools[csm]\n{e}")
        raise typer.Exit(1) from None

    lat_range = None
    if minlat is not None and maxlat is not None:
        lat_range = (minlat, maxlat)

    lon_range = None
    if minlon is not None and maxlon is not None:
        lon_range = (minlon, maxlon)

    interp_map = {
        "nearest": Interpolation.NEAREST,
        "bilinear": Interpolation.BILINEAR,
        "bicubic": Interpolation.BICUBIC,
    }
    interpolation = interp_map.get(interp.lower(), Interpolation.BICUBIC)

    _csm2map(
        input_cube=from_cube,
        output_path=to,
        map_file=map,
        projection=projection,
        resolution=resolution,
        lat_range=lat_range,
        lon_range=lon_range,
        coarse_step=step,
        dense=dense,
        validate=validate,
        clip_to_footprint=clip_to_footprint,
        shape_model=shape_model,
        spice_source=spice_source,
        interpolation=interpolation,
        profile=profile,
    )


@app.command()
def compare(
    isis_projected: Annotated[
        Path,
        typer.Argument(help="ISIS cam2map output cube.", exists=True),
    ],
    csm_projected: Annotated[
        Path,
        typer.Argument(help="csm2map output GeoTIFF.", exists=True),
    ],
):
    """Compare ISIS cam2map output with csm2map output.

    Reports pixel-level difference statistics for validation.
    """
    try:
        from isistools.csm2map.compare import compare as _compare
    except ImportError as e:
        typer.echo(f"csm2map requires extra dependencies: pip install isistools[csm]\n{e}")
        raise typer.Exit(1) from None

    try:
        _compare(isis_projected, csm_projected)
    except ValueError as e:
        typer.echo(f"[red]{e}[/red]")
        raise typer.Exit(1) from None


if __name__ == "__main__":
    app()
