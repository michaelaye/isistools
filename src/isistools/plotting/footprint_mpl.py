"""Matplotlib footprint plotting.

Provides both an interactive native window viewer and a static PNG
export suitable for reports and blog posts.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from isistools.plotting.styles import ctx_short_pid

if TYPE_CHECKING:
    import geopandas as gpd
    import pandas as pd


def _plot_footprints(
    gdf: gpd.GeoDataFrame,
    cnet_df: pd.DataFrame | None = None,
    title: str = "Footprints",
    figsize: tuple[float, float] = (12, 8),
) -> plt.Figure:
    """Create a footprint figure with optional control network overlay.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Footprints as returned by :func:`isistools.io.footprints.load_footprints`.
    cnet_df : pd.DataFrame, optional
        Control network. If provided, points are overlaid and color-coded
        by status.
    title : str
        Figure title.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, ax = plt.subplots(figsize=figsize)

    # -- Footprint polygons, colored per filename --
    filenames = gdf["filename"].unique().tolist()
    cmap = plt.colormaps["tab20"]
    colors = {fn: cmap(i % 20) for i, fn in enumerate(filenames)}

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
        legend_handles.append(Patch(edgecolor=color, facecolor=(*color[:3], 0.3),
                                    linewidth=1.5, label=ctx_short_pid(fn)))

    # -- Control point overlay --
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
            cnet_df, cube_paths=cube_paths, clock_lookup=clock_lookup,
        )
        _cnet_mpl_styles = {
            "registered": {"color": "black", "marker": "o", "markersize": 3, "alpha": 0.7},
            "unregistered": {"color": "#e74c3c", "marker": "o", "markersize": 4, "alpha": 0.7},
            "ignored": {"color": "red", "marker": "o", "markersize": 5, "alpha": 0.9},
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
        cnet_handles, _ = ax.get_legend_handles_labels()
        legend_handles.extend(cnet_handles)

    ax.legend(handles=legend_handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=11, framealpha=0.9, handlelength=2.0, handleheight=1.2,
              borderpad=0.8, labelspacing=0.6)

    xmin, ymin, xmax, ymax = gdf.total_bounds
    x_margin = (xmax - xmin) * 0.1
    ax.set_xlim(xmin - x_margin, xmax + x_margin)

    ax.set_aspect("equal")
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    fig.tight_layout()
    return fig


def footprint_png(
    gdf: gpd.GeoDataFrame,
    outpath: str | Path,
    cnet_df: pd.DataFrame | None = None,
    title: str = "Footprints",
    dpi: int = 150,
) -> Path:
    """Save a publication-ready footprint overview as PNG.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Footprints.
    outpath : path-like
        Output PNG path.
    cnet_df : pd.DataFrame, optional
        Control network overlay.
    title : str
        Figure title.
    dpi : int
        Output resolution.

    Returns
    -------
    Path
        The written file path.
    """
    import matplotlib
    matplotlib.use("Agg")

    outpath = Path(outpath)
    fig = _plot_footprints(gdf, cnet_df=cnet_df, title=title)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return outpath


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
        pass


def footprint_window(
    gdf: gpd.GeoDataFrame,
    cnet_df: pd.DataFrame | None = None,
    title: str = "Footprints",
) -> None:
    """Display footprints in a native matplotlib window with hover tooltips."""
    import matplotlib
    matplotlib.use("QtAgg")
    import mplcursors

    fig = _plot_footprints(gdf, cnet_df=cnet_df, title=title)
    fig.canvas.manager.set_window_title(title)
    _center_window(fig)

    # Add hover tooltips for interactive use
    artists_to_filename: dict = {}
    filenames = gdf["filename"].unique().tolist()
    for artist in fig.axes[0].collections:
        for fn in filenames:
            if fn not in artists_to_filename.values():
                artists_to_filename[artist] = fn
                break

    footprint_artists = list(artists_to_filename.keys())
    if footprint_artists:
        cursor = mplcursors.cursor(footprint_artists, hover=True)

        @cursor.connect("add")
        def _on_add(sel):
            sel.annotation.set_text(artists_to_filename.get(sel.artist, ""))
            sel.annotation.get_bbox_patch().set(fc="white", alpha=0.9)

    plt.show()
