"""Command-line interface for isistools.

All commands launch Panel apps in the browser.

Usage::

    isistools mosaic cubes.lis --cnet control.net
    isistools tiepoints cubes.lis control.net
    isistools footprints cubes.lis
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

app = typer.Typer(
    name="isistools",
    help="Python review tools for ISIS3 coregistration workflows.",
    no_args_is_help=True,
)


@app.command()
def mosaic(
    cubelist: Path = typer.Argument(
        ..., help="Cube list file (one cube path per line)", exists=True,
    ),
    cnet: Optional[Path] = typer.Option(
        None, "--cnet", "-c", help="Control network file (.net)",
    ),
    port: int = typer.Option(0, "--port", "-p", help="Server port (0=auto)"),
    no_browser: bool = typer.Option(False, "--no-browser", help="Don't open browser"),
):
    """Review mosaic footprints and image content (Qmos replacement)."""
    from isistools.apps.mosaic_review import MosaicReview

    typer.echo(f"Loading cubes from {cubelist}...")
    review = MosaicReview(cube_list=cubelist, cnet_path=cnet)
    typer.echo("Starting server...")
    review.serve(port=port, show=not no_browser)


@app.command()
def tiepoints(
    cubelist: Path = typer.Argument(
        ..., help="Cube list file", exists=True,
    ),
    cnet: Path = typer.Argument(
        ..., help="Control network file (.net)", exists=True,
    ),
    port: int = typer.Option(0, "--port", "-p", help="Server port (0=auto)"),
    no_browser: bool = typer.Option(False, "--no-browser", help="Don't open browser"),
):
    """Review tie points between image pairs (Qnet replacement)."""
    from isistools.apps.tiepoint_review import TiepointReview

    typer.echo(f"Loading {cubelist} with {cnet}...")
    review = TiepointReview(cube_list=cubelist, cnet_path=cnet)
    typer.echo("Starting server...")
    review.serve(port=port, show=not no_browser)


@app.command()
def footprints(
    cubelist: Path = typer.Argument(
        ..., help="Cube list file", exists=True,
    ),
    cnet: Optional[Path] = typer.Option(
        None, "--cnet", "-c", help="Optional control network overlay",
    ),
    port: int = typer.Option(0, "--port", "-p", help="Server port (0=auto)"),
    no_browser: bool = typer.Option(False, "--no-browser", help="Don't open browser"),
):
    """Quick footprint map viewer."""
    import panel as pn

    from isistools.io.footprints import load_footprints
    from isistools.plotting.footprint_map import footprint_map, footprint_map_with_cnet

    typer.echo(f"Loading footprints from {cubelist}...")
    gdf = load_footprints(cubelist)
    typer.echo(f"Loaded {len(gdf)} footprints")

    if cnet is not None:
        from isistools.io.controlnet import load_cnet
        from isistools.plotting.cnet_overlay import cnet_to_geodataframe

        cnet_df = load_cnet(cnet)
        cnet_gdf = cnet_to_geodataframe(cnet_df)
        plot = footprint_map_with_cnet(gdf, cnet_gdf)
    else:
        plot = footprint_map(gdf)

    pn.serve(pn.pane.HoloViews(plot), port=port, show=not no_browser, title="Footprints")


@app.command()
def cnet_info(
    cnet: Path = typer.Argument(
        ..., help="Control network file (.net)", exists=True,
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


if __name__ == "__main__":
    app()
