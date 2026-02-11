"""Control network point overlay visualization.

Renders control network tie points on top of images (in sample/line
space) or maps (in lon/lat space), using the improved color scheme
defined in styles.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import geopandas as gpd
import holoviews as hv
import hvplot.pandas  # noqa: F401
import numpy as np
import pandas as pd
from shapely.geometry import Point

from isistools.plotting.styles import CNET_POINT_STYLES, STATUS_COLOR_MAP

if TYPE_CHECKING:
    pass

hv.extension("bokeh")


def cnet_points_image(
    cnet_df: pd.DataFrame,
    serial_number: str | None = None,
    hover_cols: list[str] | None = None,
) -> hv.Element:
    """Overlay control points on an image in sample/line space.

    Parameters
    ----------
    cnet_df : pd.DataFrame
        Control network (one row per measure).
    serial_number : str, optional
        Filter to only show measures for this image.
    hover_cols : list of str, optional
        Columns to show on hover. Defaults to pointId, status,
        residual_magnitude.

    Returns
    -------
    holoviews.Element
        Scatter overlay in sample/line coordinates.
    """
    df = cnet_df.copy()

    if serial_number is not None:
        df = df[df["serialnumber"] == serial_number]

    if df.empty:
        return hv.Points([])

    if hover_cols is None:
        hover_cols = [c for c in ["pointId", "status", "residual_magnitude",
                                   "measureType", "pointType"]
                      if c in df.columns]

    # Build per-status overlays for proper styling
    overlays = []
    for status, style in CNET_POINT_STYLES.items():
        if status == "selected":
            continue  # interactive-only
        subset = df[df["status"] == status]
        if subset.empty:
            continue

        scatter = subset.hvplot.scatter(
            x="sample",
            y="line",
            hover_cols=hover_cols,
            color=style["color"],
            alpha=style["alpha"],
            size=style["size"] * 10,  # hvplot size is area-based
            label=f"{status} ({len(subset)})",
        )
        overlays.append(scatter)

    if not overlays:
        return hv.Points([])

    result = overlays[0]
    for overlay in overlays[1:]:
        result = result * overlay

    return result


def cnet_points_map(
    cnet_gdf: gpd.GeoDataFrame,
    hover_cols: list[str] | None = None,
) -> hv.Element:
    """Overlay control points on a map in lon/lat space.

    Parameters
    ----------
    cnet_gdf : gpd.GeoDataFrame
        Control points with lon/lat geometry and ``status`` column.
    hover_cols : list of str, optional
        Columns to show on hover.

    Returns
    -------
    holoviews.Element
        Points overlay for map plots.
    """
    if hover_cols is None:
        hover_cols = [c for c in ["pointId", "status", "residual_magnitude",
                                   "n_measures"]
                      if c in cnet_gdf.columns]

    overlays = []
    for status, style in CNET_POINT_STYLES.items():
        if status == "selected":
            continue
        subset = cnet_gdf[cnet_gdf["status"] == status]
        if subset.empty:
            continue

        scatter = subset.hvplot.points(
            geo=True,
            hover_cols=hover_cols,
            color=style["color"],
            alpha=style["alpha"],
            size=style["size"] * 10,
            label=f"{status} ({len(subset)})",
        )
        overlays.append(scatter)

    if not overlays:
        return hv.Points([])

    result = overlays[0]
    for overlay in overlays[1:]:
        result = result * overlay

    return result


def cnet_to_geodataframe(cnet_df: pd.DataFrame) -> gpd.GeoDataFrame:
    """Convert a control network DataFrame to a GeoDataFrame.

    Groups measures by point, takes the average adjusted lat/lon
    (or apriori if adjusted is not available), and creates point
    geometries.

    Parameters
    ----------
    cnet_df : pd.DataFrame
        Control network with one row per measure.

    Returns
    -------
    gpd.GeoDataFrame
        One row per control point with lon/lat geometry.
    """
    # Group by point to get per-point info
    lon_col = None
    lat_col = None

    # Try adjusted coordinates first, fall back to apriori
    for lon_candidate in ["adjustedX", "adjustedLon", "aprioriX", "aprioriLon"]:
        if lon_candidate in cnet_df.columns:
            lon_col = lon_candidate
            break
    for lat_candidate in ["adjustedY", "adjustedLat", "aprioriY", "aprioriLat"]:
        if lat_candidate in cnet_df.columns:
            lat_col = lat_candidate
            break

    if lon_col is None or lat_col is None:
        raise ValueError(
            "Cannot find longitude/latitude columns in control network. "
            f"Available columns: {list(cnet_df.columns)}"
        )

    # Aggregate per point
    point_groups = cnet_df.groupby("pointId")

    records = []
    for point_id, group in point_groups:
        lon = group[lon_col].iloc[0]
        lat = group[lat_col].iloc[0]

        # Determine overall point status
        statuses = group["status"].unique()
        if "ignored" in statuses and len(statuses) == 1:
            status = "ignored"
        elif "registered" in statuses:
            status = "registered"
        else:
            status = "unregistered"

        records.append({
            "pointId": point_id,
            "geometry": Point(float(lon), float(lat)),
            "status": status,
            "n_measures": len(group),
            "residual_magnitude": group["residual_magnitude"].mean(),
            "pointType": group.get("pointType", pd.Series("Unknown")).iloc[0],
        })

    gdf = gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")
    return gdf


def cnet_residual_vectors(
    cnet_df: pd.DataFrame,
    serial_number: str | None = None,
    scale: float = 10.0,
) -> hv.Element:
    """Plot residual vectors as arrows on an image.

    Parameters
    ----------
    cnet_df : pd.DataFrame
        Control network measures.
    serial_number : str, optional
        Filter to specific image.
    scale : float
        Scaling factor for residual vectors (pixels per residual unit).

    Returns
    -------
    holoviews.Element
        Vectorfield or segments overlay.
    """
    df = cnet_df.copy()
    if serial_number is not None:
        df = df[df["serialnumber"] == serial_number]

    df = df[df["status"] == "registered"].copy()
    if df.empty:
        return hv.Segments([])

    # Create segments from (sample, line) to (sample + res*scale, line + res*scale)
    df["sample_end"] = df["sample"] + df.get("residualSample", 0.0) * scale
    df["line_end"] = df["line"] + df.get("residualLine", 0.0) * scale

    segments = hv.Segments(
        df, kdims=["sample", "line", "sample_end", "line_end"]
    ).opts(
        color="#e74c3c",
        line_width=1.5,
        alpha=0.8,
    )

    return segments
