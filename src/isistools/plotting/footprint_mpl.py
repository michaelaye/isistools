"""Native matplotlib footprint viewer.

Lightweight alternative to the browser-based hvplot viewer. Launches
a native window with zoom/pan toolbar and hover tooltips showing the
CTX product ID (filename).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import matplotlib
matplotlib.use("QtAgg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import mplcursors

from isistools.plotting.styles import ctx_short_pid

if TYPE_CHECKING:
    import geopandas as gpd
    import pandas as pd


def _center_window(fig) -> None:
    """Attempt to center the matplotlib window on screen."""
    try:
        backend = plt.get_backend()
        manager = fig.canvas.manager
        if "Tk" in backend:
            window = manager.window
            window.update_idletasks()
            sw = window.winfo_screenwidth()
            sh = window.winfo_screenheight()
            ww = window.winfo_width()
            wh = window.winfo_height()
            window.geometry(f"+{(sw - ww) // 2}+{(sh - wh) // 2}")
        elif "Qt" in backend:
            window = manager.window
            screen = window.screen().geometry()
            size = window.frameGeometry()
            window.move(
                (screen.width() - size.width()) // 2,
                (screen.height() - size.height()) // 2,
            )
    except Exception:
        pass  # Not all backends support window positioning


def footprint_window(
    gdf: gpd.GeoDataFrame,
    cnet_df: pd.DataFrame | None = None,
    title: str = "Footprints",
) -> None:
    """Display footprints in a native matplotlib window.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Footprints as returned by :func:`isistools.io.footprints.load_footprints`.
    cnet_df : pd.DataFrame, optional
        Control network. If provided, points are overlaid and color-coded
        by status. Uses ISIS campt for precise coordinate conversion on
        pre-jigsaw networks without ground coordinates.
    title : str
        Window title.
    """
    fig, ax = plt.subplots(figsize=(12, 8))
    fig.canvas.manager.set_window_title(title)
    _center_window(fig)

    # -- Footprint polygons, colored per filename --
    filenames = gdf["filename"].unique().tolist()
    cmap = plt.colormaps["tab20"]
    colors = {fn: cmap(i % 20) for i, fn in enumerate(filenames)}

    # Plot each filename separately so mplcursors can identify them
    artists_to_filename: dict = {}
    legend_handles = []
    for fn in filenames:
        sub = gdf[gdf["filename"] == fn]
        color = colors[fn]
        sub.plot(
            ax=ax,
            facecolor=(*color[:3], 0.3),
            edgecolor=color,
            linewidth=1.5,
        )
        # gdf.plot returns the axes; the new artists are the last PathCollection(s)
        for artist in ax.collections:
            if artist not in artists_to_filename:
                artists_to_filename[artist] = fn
        # Legend entry with short product ID (first 18 chars)
        legend_handles.append(Patch(edgecolor=color, facecolor=(*color[:3], 0.3),
                                    linewidth=1.5, label=ctx_short_pid(fn)))

    # Hover tooltip for footprints — shows filename on hover
    footprint_artists = list(artists_to_filename.keys())
    if footprint_artists:
        cursor = mplcursors.cursor(footprint_artists, hover=True)

        @cursor.connect("add")
        def _on_add(sel):
            sel.annotation.set_text(artists_to_filename.get(sel.artist, ""))
            sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)

    # -- Control point overlay --
    if cnet_df is not None:
        from pathlib import Path

        from isistools.plotting.cnet_overlay import cnet_to_geodataframe

        cube_paths = gdf["path"].tolist() if "path" in gdf.columns else None
        # Reuse clock counts already extracted during load_footprints
        clock_lookup = None
        if "clock" in gdf.columns and "path" in gdf.columns:
            clock_lookup = {
                row["clock"]: Path(row["path"])
                for _, row in gdf[["clock", "path"]].iterrows()
                if row["clock"]
            }
        cnet_gdf = cnet_to_geodataframe(
            cnet_df, cube_paths=cube_paths, clock_lookup=clock_lookup,
        )
        _cnet_mpl_styles = {
            "registered": {"color": "black", "marker": "x", "markersize": 16, "alpha": 0.9},
            "unregistered": {"color": "#e74c3c", "marker": "o", "markersize": 4, "alpha": 0.7},
            "ignored": {"color": "#95a5a6", "marker": "o", "markersize": 3, "alpha": 0.4},
        }
        for status, style in _cnet_mpl_styles.items():
            sub = cnet_gdf[cnet_gdf["status"] == status]
            if sub.empty:
                continue
            sub.plot(
                ax=ax,
                **style,
                label=f"{status} ({len(sub)})",
            )
        # Combine cnet status entries with footprint legend
        cnet_handles, cnet_labels = ax.get_legend_handles_labels()
        legend_handles.extend(cnet_handles)

    ax.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=11, framealpha=0.9, handlelength=2.0, handleheight=1.2,
              borderpad=0.8, labelspacing=0.6)

    # Constrain x-axis to data range +/- 10%
    xmin, ymin, xmax, ymax = gdf.total_bounds
    x_margin = (xmax - xmin) * 0.1
    ax.set_xlim(xmin - x_margin, xmax + x_margin)

    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    plt.tight_layout()
    plt.show()
