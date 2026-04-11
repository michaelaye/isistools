"""Command-line interface for isistools.

Usage::

    isistools mosaic cubes.lis --cnet control.net
    isistools tiepoints cubes.lis control.net
    isistools footprints cubes.lis
    isistools csm2map input.cub output.tif --map equi.map
"""

import os
from pathlib import Path
from typing import Annotated, Optional

import typer


def _ensure_cwd() -> None:
    """Ensure the current working directory exists.

    The ``param`` library calls ``os.getcwd()`` at import time (during class
    definition of ``resolve_path``).  On macOS the CWD can transiently vanish
    (e.g. when launched from a deleted tmpdir), causing an immediate
    ``FileNotFoundError`` on ``import panel`` / ``import holoviews``.

    Call this before any panel/holoviews/param import.
    """
    try:
        os.getcwd()
    except FileNotFoundError:
        os.chdir(Path.home())


app = typer.Typer(
    name="isistools",
    help="Python review tools for ISIS3 coregistration workflows.",
    no_args_is_help=True,
)


@app.command()
def mosaic(
    cubelist: Path = typer.Argument(
        ...,
        help="Cube list file (one cube path per line)",
        exists=True,
    ),
    cnet: Optional[Path] = typer.Option(
        None,
        "--cnet",
        "-c",
        help="Control network file (.net)",
    ),
    port: int = typer.Option(0, "--port", "-p", help="Server port (0=auto)"),
    no_browser: bool = typer.Option(False, "--no-browser", help="Don't open browser"),
):
    """Review mosaic footprints and image content (Qmos replacement)."""
    _ensure_cwd()
    from isistools.apps.mosaic_review import MosaicReview

    typer.echo(f"Loading cubes from {cubelist}...")
    review = MosaicReview(cube_list=cubelist, cnet_path=cnet)
    typer.echo("Starting server...")
    review.serve(port=port, show=not no_browser)


@app.command()
def tiepoints(
    cubelist: Path = typer.Argument(
        ...,
        help="Cube list file",
        exists=True,
    ),
    cnet: Path = typer.Argument(
        ...,
        help="Control network file (.net)",
        exists=True,
    ),
    port: int = typer.Option(0, "--port", "-p", help="Server port (0=auto)"),
    no_browser: bool = typer.Option(False, "--no-browser", help="Don't open browser"),
):
    """Review tie points between image pairs (Qnet replacement)."""
    _ensure_cwd()
    from isistools.apps.tiepoint_review import TiepointReview

    typer.echo(f"Loading {cubelist} with {cnet}...")
    review = TiepointReview(cube_list=cubelist, cnet_path=cnet)
    typer.echo("Starting server...")
    review.serve(port=port, show=not no_browser)


@app.command()
def footprints(
    cubelist: Path = typer.Argument(
        ...,
        help="Cube list file",
        exists=True,
    ),
    cnet: Optional[Path] = typer.Option(
        None,
        "--cnet",
        "-c",
        help="Optional control network overlay",
    ),
    png: bool = typer.Option(
        False,
        "--png",
        help="Save static PNG instead of launching viewer",
    ),
    png_path: Optional[Path] = typer.Option(
        None,
        "--png-path",
        help="PNG output path (default: footprints_overview.png)",
    ),
    dpi: int = typer.Option(150, "--dpi", help="PNG resolution (only with --png)"),
    title: Optional[str] = typer.Option(
        None,
        "--title",
        "-t",
        help="Figure title (default: cubelist filename)",
    ),
    win: bool = typer.Option(
        False,
        "--win",
        help="Native matplotlib window instead of browser",
    ),
    port: int = typer.Option(0, "--port", "-p", help="Server port (0=auto)"),
    no_browser: bool = typer.Option(False, "--no-browser", help="Don't open browser"),
):
    """Quick footprint map viewer, or static PNG export with --png."""
    from isistools.io.footprints import load_footprints

    typer.echo(f"Loading footprints from {cubelist}...")
    gdf = load_footprints(cubelist)
    typer.echo(f"Loaded {len(gdf)} footprints")

    fig_title = title or cubelist.stem

    cnet_df = None
    if cnet is not None:
        from isistools.io.controlnet import load_cnet

        cnet_df = load_cnet(cnet)

    if png or png_path:
        from isistools.plotting.footprint_mpl import footprint_png

        outfile = png_path or Path("footprints_overview.png")
        out = footprint_png(gdf, outfile, cnet_df=cnet_df, title=fig_title, dpi=dpi)
        typer.echo(f"Saved: {out}")
    elif win:
        from isistools.plotting.footprint_mpl import footprint_window

        footprint_window(gdf, cnet_df=cnet_df, title=fig_title)
    else:
        _ensure_cwd()
        import panel as pn

        from isistools.plotting.footprint_map import footprint_map, footprint_map_with_cnet

        if cnet_df is not None:
            from isistools.plotting.cnet_overlay import cnet_to_geodataframe

            cube_paths = gdf["path"].tolist() if "path" in gdf.columns else None
            clock_lookup = None
            if "clock" in gdf.columns and "path" in gdf.columns:
                clock_lookup = {
                    row["clock"]: Path(row["path"])
                    for _, row in gdf[["clock", "path"]].iterrows()
                    if row["clock"]
                }
            cnet_gdf = cnet_to_geodataframe(
                cnet_df,
                cube_paths=cube_paths,
                clock_lookup=clock_lookup,
            )
            plot = footprint_map_with_cnet(gdf, cnet_gdf)
        else:
            plot = footprint_map(gdf)

        pn.serve(
            pn.pane.HoloViews(plot),
            port=port,
            show=not no_browser,
            title="Footprints",
        )


@app.command()
def footprintinit(
    cubelist: Path = typer.Argument(
        ...,
        help="Cube list file (one cube path per line)",
        exists=True,
    ),
    jobs: int = typer.Option(
        4,
        "--jobs",
        "-j",
        help="Number of parallel workers",
    ),
):
    """Run ISIS footprintinit on all cubes in a list file."""
    import subprocess
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from isistools.io.footprints import read_cube_list

    cubes = read_cube_list(cubelist)
    typer.echo(f"Running footprintinit on {len(cubes)} cubes ({jobs} workers)...")

    def _run(cube: Path) -> tuple[Path, bool, str]:
        result = subprocess.run(
            ["footprintinit", f"from={cube}"],
            capture_output=True,
            text=True,
        )
        ok = result.returncode == 0
        msg = result.stderr.strip() if not ok else ""
        return cube, ok, msg

    failed = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(_run, c): c for c in cubes}
        for future in as_completed(futures):
            cube, ok, msg = future.result()
            if ok:
                typer.echo(f"  OK: {cube.name}")
            else:
                typer.echo(f"  FAIL: {cube.name} — {msg}")
                failed.append(cube)

    if failed:
        typer.echo(f"\n{len(failed)}/{len(cubes)} failed.")
        raise typer.Exit(1)
    typer.echo(f"\nAll {len(cubes)} cubes done.")


@app.command()
def spiceinit(
    cubelist: Path = typer.Argument(
        ...,
        help="Cube list file (one cube path per line)",
        exists=True,
    ),
    web: bool = typer.Option(
        True,
        "--web/-W",
        "-w",
        help="Use web=yes for SPICE kernel retrieval (default: on)",
    ),
    jobs: int = typer.Option(
        4,
        "--jobs",
        "-j",
        help="Number of parallel workers",
    ),
):
    """Run ISIS spiceinit on all cubes in a list file."""
    import subprocess
    from concurrent.futures import ThreadPoolExecutor, as_completed

    from isistools.io.footprints import read_cube_list

    cubes = read_cube_list(cubelist)
    typer.echo(f"Running spiceinit on {len(cubes)} cubes ({jobs} workers)...")

    def _run(cube: Path) -> tuple[Path, bool, str]:
        cmd = ["spiceinit", f"from={cube}"]
        if web:
            cmd.append("web=yes")
        result = subprocess.run(cmd, capture_output=True, text=True)
        ok = result.returncode == 0
        msg = result.stderr.strip() if not ok else ""
        return cube, ok, msg

    failed = []
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = {pool.submit(_run, c): c for c in cubes}
        for future in as_completed(futures):
            cube, ok, msg = future.result()
            if ok:
                typer.echo(f"  OK: {cube.name}")
            else:
                typer.echo(f"  FAIL: {cube.name} — {msg}")
                failed.append(cube)

    if failed:
        typer.echo(f"\n{len(failed)}/{len(cubes)} failed.")
        raise typer.Exit(1)
    typer.echo(f"\nAll {len(cubes)} cubes done.")


@app.command()
def overlaps(
    cubelist: Path = typer.Argument(
        ...,
        help="Cube list file (one cube path per line)",
        exists=True,
    ),
    output: Optional[Path] = typer.Option(
        None,
        "--output",
        "-o",
        help="Output overlap list path (default: <cubelist_dir>/overlap_list.lis)",
    ),
    png: bool = typer.Option(
        False,
        "--png",
        help="Save a PNG plot of the overlap polygons",
    ),
    png_path: Optional[Path] = typer.Option(
        None,
        "--png-path",
        help="PNG output path (default: overlaps.png)",
    ),
    dpi: int = typer.Option(150, "--dpi", help="PNG resolution"),
    gpkg: Optional[Path] = typer.Option(
        None,
        "--gpkg",
        help="Export overlap polygons to GeoPackage",
    ),
):
    """Run findimageoverlaps and extract overlap polygons as GeoDataFrame.

    Runs ISIS findimageoverlaps, then parses the output into a GeoDataFrame
    with polygon geometries for each overlap zone. Optionally exports to
    GeoPackage or PNG.
    """
    import subprocess

    from isistools.io.overlaps import parse_overlap_list

    overlap_out = output or (cubelist.parent / "overlap_list.lis")

    typer.echo(f"Running findimageoverlaps on {cubelist}...")
    result = subprocess.run(
        [
            "findimageoverlaps",
            f"fromlist={cubelist}",
            f"overlaplist={overlap_out}",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        typer.echo(f"ERROR: findimageoverlaps failed:\n{result.stderr}")
        raise typer.Exit(1)

    typer.echo(f"Overlap list written to: {overlap_out}")
    typer.echo("Parsing overlap polygons...")

    gdf = parse_overlap_list(overlap_out)

    # Summary table
    typer.echo(f"\n{'Zone':<50s} {'Type':<16s} {'Area (deg²)':<12s}")
    typer.echo("-" * 78)
    for _, row in gdf.iterrows():
        serials = row["serials"]
        # Truncate long serial strings
        if len(serials) > 48:
            serials = serials[:45] + "..."
        typer.echo(f"{serials:<50s} {row['zone_type']:<16s} {row['area_deg2']:<12.6f}")

    overlap_only = gdf[gdf["n_images"] >= 2]
    typer.echo(
        f"\n{len(overlap_only)} overlap zones, "
        f"{len(gdf) - len(overlap_only)} individual footprints"
    )

    if gpkg:
        # Convert list column to string for GeoPackage compatibility
        export = gdf.copy()
        export["images"] = export["images"].apply(lambda x: ",".join(x))
        export.to_file(gpkg, driver="GPKG")
        typer.echo(f"Saved GeoPackage: {gpkg}")

    if png or png_path:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch

        colors = {
            "footprint": "#cccccc",
            "2-way overlap": "#4a90d9",
        }
        # Catch N-way overlaps
        default_overlap_color = "#e74c3c"

        fig, ax = plt.subplots(1, 1, figsize=(10, 8))

        for _, row in gdf.iterrows():
            zt = row["zone_type"]
            color = colors.get(zt, default_overlap_color)
            alpha = 0.3 if zt == "footprint" else 0.5
            gpd.GeoDataFrame([row], geometry="geometry", crs="EPSG:4326").plot(
                ax=ax,
                color=color,
                alpha=alpha,
                edgecolor="black",
                linewidth=0.8,
            )
            centroid = row.geometry.centroid
            ax.annotate(
                row["serials"].replace("MRO/CTX/", "").replace(",", "\n"),
                xy=(centroid.x, centroid.y),
                ha="center",
                va="center",
                fontsize=6,
            )

        handles = []
        seen = set()
        for _, row in gdf.iterrows():
            zt = row["zone_type"]
            if zt not in seen:
                color = colors.get(zt, default_overlap_color)
                handles.append(Patch(facecolor=color, alpha=0.5, edgecolor="black", label=zt))
                seen.add(zt)
        ax.legend(handles=handles, loc="upper left")
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
        ax.set_title(f"Overlap Zones — {cubelist.name}")
        ax.set_aspect("equal")

        out = png_path or Path("overlaps.png")
        fig.savefig(out, dpi=dpi, bbox_inches="tight")
        typer.echo(f"Saved PNG: {out}")
        plt.close(fig)


@app.command()
def cnet_info(
    cnet: Path = typer.Argument(
        ...,
        help="Control network file (.net)",
        exists=True,
    ),
):
    """Print control network summary statistics."""
    from isistools.io.controlnet import cnet_summary, load_cnet

    df = load_cnet(cnet)
    stats = cnet_summary(df)

    typer.echo(f"Control Network: {cnet.name}")
    typer.echo(f"  Points:        {stats['n_points']}")
    typer.echo(f"  Measures:      {stats['n_measures']}")
    typer.echo(f"  Images:        {stats['n_images']}")
    typer.echo(f"  Registered:    {stats['n_registered']}")
    typer.echo(f"  Unregistered:  {stats['n_unregistered']}")
    typer.echo(f"  Ignored:       {stats['n_ignored']}")
    typer.echo(f"  Mean Residual: {stats['mean_residual']:.4f}")
    typer.echo(f"  Max Residual:  {stats['max_residual']:.4f}")


@app.command()
def csm2map(
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
        typer.Option("--resolution", "-r", help="Pixel resolution in meters/pixel."),
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
        typer.Option("--step", "-s", help="Coarse grid step (pixels). Smaller=more accurate."),
    ] = 16,
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
            help="Clip output to the footprint polygon stored by footprintinit "
            "(ISIS cam2map compatibility mode). Default: full camera coverage.",
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
    interp: Annotated[
        str,
        typer.Option("--interp", "-i", help="Interpolation: nearest, bilinear, bicubic."),
    ] = "bicubic",
):
    """Map-project an ISIS cube using a CSM camera model (ISIS cam2map replacement).

    By default produces MORE coverage than ISIS cam2map because it uses the
    true camera model instead of the coarse footprint polygon from
    footprintinit. Use --clip-to-footprint for exact ISIS compatibility.

    By default the shape model is read from the cube label, matching ISIS
    cam2map. Use --shape-model ellipsoid to disable DEM lookups (faster but
    less accurate over topography), or --shape-model PATH to use a custom DEM.

    \b
    Examples:
      isistools csm2map input.cub output.tif --map equi.map
      isistools csm2map input.cub output.tif -r 6.0
      isistools csm2map input.cub output.tif --map equi.map --dense --validate
      isistools csm2map input.cub output.tif --map equi.map --clip-to-footprint
      isistools csm2map input.cub output.tif --map equi.map --shape-model ellipsoid
    """
    try:
        from isistools.processing.project import project
        from isistools.processing.resample import Interpolation
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

    project(
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
        interpolation=interpolation,
    )


@app.command()
def csm2map_compare(
    isis_projected: Annotated[
        Path,
        typer.Argument(help="ISIS cam2map output cube.", exists=True),
    ],
    csm_projected: Annotated[
        Path,
        typer.Argument(help="isistools csm2map output GeoTIFF.", exists=True),
    ],
):
    """Compare ISIS cam2map output with isistools csm2map output.

    Reports pixel-level difference statistics for validation.
    """
    import numpy as np
    import rasterio
    from rich.console import Console

    from isistools.io.cubes import read_isis_cube_raw

    console = Console()

    console.print("[bold]Loading ISIS projected cube[/bold]")
    isis_data, _ = read_isis_cube_raw(isis_projected)

    console.print("[bold]Loading CSM projected GeoTIFF[/bold]")
    with rasterio.open(str(csm_projected)) as src:
        csm_data = src.read(1).astype(np.float32)

    # Check shapes
    if isis_data.shape != csm_data.shape:
        console.print(f"[red]Shape mismatch: ISIS {isis_data.shape} vs CSM {csm_data.shape}[/red]")
        console.print("Comparison requires matching grid parameters.")
        raise typer.Exit(1)

    # Compare only where both have valid data
    isis_valid = np.isfinite(isis_data) & (isis_data != 0)
    csm_valid = np.isfinite(csm_data) & (csm_data != 0)
    both_valid = isis_valid & csm_valid

    n_both = int(np.sum(both_valid))
    n_isis_only = int(np.sum(isis_valid & ~csm_valid))
    n_csm_only = int(np.sum(csm_valid & ~isis_valid))

    console.print(f"\n  Both valid: {n_both:,}")
    console.print(f"  ISIS-only:  {n_isis_only:,}")
    console.print(f"  CSM-only:   {n_csm_only:,}")

    if n_both == 0:
        console.print("[red]No overlapping valid pixels![/red]")
        raise typer.Exit(1)

    diff = csm_data[both_valid] - isis_data[both_valid]
    console.print("\n  [bold]Difference statistics (CSM - ISIS):[/bold]")
    console.print(f"  Mean:   {np.mean(diff):.4f}")
    console.print(f"  Median: {np.median(diff):.4f}")
    console.print(f"  Std:    {np.std(diff):.4f}")
    console.print(f"  Min:    {np.min(diff):.4f}")
    console.print(f"  Max:    {np.max(diff):.4f}")

    # Percentage of pixels within DN thresholds
    for threshold in [0.01, 0.1, 1.0, 5.0]:
        pct = 100 * np.sum(np.abs(diff) < threshold) / n_both
        console.print(f"  |diff| < {threshold}: {pct:.1f}%")


if __name__ == "__main__":
    app()
