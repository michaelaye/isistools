"""Interactive footprint map visualization.

Displays ISIS cube footprints on an interactive map using hvplot,
with hover information and click-to-select functionality.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import holoviews as hv
import hvplot.pandas  # noqa: F401

from isistools.plotting.styles import ctx_short_pid

if TYPE_CHECKING:
    import geopandas as gpd

hv.extension("bokeh")


def footprint_map(
    gdf: gpd.GeoDataFrame,
    hover_cols: list[str] | None = None,
    title: str = "Image Footprints",
    width: int = 1200,
    height: int = 800,
) -> hv.Element:
    """Create an interactive footprint overview map.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Footprints as returned by :func:`isistools.io.footprints.load_footprints`.
    hover_cols : list of str, optional
        Columns to show on hover. Defaults to filename, instrument, start_time.
    title : str
        Plot title.
    width, height : int
        Plot dimensions in pixels.

    Returns
    -------
    holoviews.Element
        Interactive map plot (renders in notebooks and Panel apps).
    """
    if hover_cols is None:
        hover_cols = [c for c in ["start_time"]
                      if c in gdf.columns]

    gdf = gdf.copy()
    gdf["short_pid"] = gdf["filename"].map(ctx_short_pid)
    if "start_time" in gdf.columns:
        gdf["start_time"] = gdf["start_time"].astype(str).str[:19]

    def _plot_tweaks(plot, element):
        from bokeh.models import BoxZoomTool, Legend
        for legend in plot.state.select({"type": Legend}):
            legend.spacing = 20
            legend.padding = 15
            legend.label_standoff = 8
        for tool in plot.state.select({"type": BoxZoomTool}):
            tool.match_aspect = False

    plot = gdf.hvplot(
        c="short_pid",
        hover_cols=hover_cols,
        fill_alpha=0.3,
        line_width=2,
        legend="bottom",
        legend_cols=3,
        width=width,
        height=height,
        title=title,
        tools=["tap", "hover", "wheel_zoom", "pan", "reset"],
        xlabel="Longitude",
        ylabel="Latitude",
        fontsize={
            "xlabel": "20pt",
            "ylabel": "20pt",
            "xticks": "14pt",
            "yticks": "14pt",
            "legend": "15pt",
        },
    ).opts(hooks=[_plot_tweaks], data_aspect=1)

    return plot


def footprint_map_with_cnet(
    gdf: gpd.GeoDataFrame,
    cnet_gdf: gpd.GeoDataFrame,
    **kwargs,
) -> hv.Element:
    """Footprint map with control network points overlaid.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Image footprints.
    cnet_gdf : gpd.GeoDataFrame
        Control points as a GeoDataFrame with lon/lat geometry
        and a ``status`` column.
    **kwargs
        Passed to :func:`footprint_map`.

    Returns
    -------
    holoviews.Element
        Combined footprint + points overlay.
    """
    from isistools.plotting.cnet_overlay import cnet_points_map

    base = footprint_map(gdf, **kwargs)
    points = cnet_points_map(cnet_gdf)

    return base * points
