"""Interactive footprint map visualization.

Displays ISIS cube footprints on an interactive map using hvplot,
with hover information and click-to-select functionality.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import holoviews as hv
import hvplot.pandas  # noqa: F401

from isistools.plotting.styles import FOOTPRINT_STYLES

if TYPE_CHECKING:
    import geopandas as gpd

hv.extension("bokeh")


def footprint_map(
    gdf: gpd.GeoDataFrame,
    hover_cols: list[str] | None = None,
    title: str = "Image Footprints",
    width: int = 800,
    height: int = 500,
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
        hover_cols = [c for c in ["filename", "instrument", "start_time", "target"]
                      if c in gdf.columns]

    from bokeh.palettes import Category20_20

    filenames = gdf["filename"].unique().tolist()
    colors = {fn: Category20_20[i % 20] for i, fn in enumerate(filenames)}

    overlays = []
    for fn in filenames:
        sub = gdf[gdf["filename"] == fn]
        p = sub.hvplot(
            geo=True,
            hover_cols=hover_cols,
            fill_color=colors[fn],
            fill_alpha=0.3,
            line_color=colors[fn],
            line_width=2,
            label=fn,
        )
        overlays.append(p)

    plot = hv.Overlay(overlays).opts(
        hv.opts.Polygons(muted_alpha=0.05, muted_fill_alpha=0.02),
        hv.opts.Overlay(
            width=width, height=height, title=title,
            tools=["tap", "hover", "wheel_zoom", "pan", "reset"],
            legend_position="bottom",
            legend_opts={
                "click_policy": "mute",
                "orientation": "vertical",
                "ncols": 2,
                "label_text_font_size": "14px",
                "spacing": 2,
                "padding": 5,
            },
        ),
    )

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
