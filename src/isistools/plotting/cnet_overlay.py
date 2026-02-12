"""Control network point overlay visualization.

Renders control network tie points on top of images (in sample/line
space) or maps (in lon/lat space), using the improved color scheme
defined in styles.py.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point

from isistools.plotting.styles import CNET_POINT_STYLES, STATUS_COLOR_MAP

if TYPE_CHECKING:
    import holoviews as hv


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
    import holoviews as hv
    import hvplot.pandas  # noqa: F401

    hv.extension("bokeh")

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
    import holoviews as hv
    import hvplot.pandas  # noqa: F401

    hv.extension("bokeh")

    if hover_cols is None:
        hover_cols = [c for c in ["status", "residual_magnitude",
                                   "n_measures"]
                      if c in cnet_gdf.columns]

    # Styles matching the --win matplotlib path
    _map_styles = {
        "registered": {"color": "black", "marker": "x", "size": 160, "alpha": 0.9},
        "unregistered": {"color": "#e74c3c", "marker": "circle", "size": 40, "alpha": 0.7},
        "ignored": {"color": "#95a5a6", "marker": "circle", "size": 30, "alpha": 0.4},
    }

    overlays = []
    for status, style in _map_styles.items():
        subset = cnet_gdf[cnet_gdf["status"] == status]
        if subset.empty:
            continue

        scatter = subset[["geometry"]].hvplot.points(
            hover=False,
            color=style["color"],
            marker=style["marker"],
            alpha=style["alpha"],
            size=style["size"],
            label=f"{status} ({len(subset)})",
        )
        overlays.append(scatter)

    if not overlays:
        return hv.Points([])

    result = overlays[0]
    for overlay in overlays[1:]:
        result = result * overlay

    return result


def _has_ground_coords(cnet_df: pd.DataFrame) -> bool:
    """Check if the control network has non-zero ground coordinates."""
    for col in ["adjustedX", "adjustedLon", "aprioriX", "aprioriLon"]:
        if col in cnet_df.columns and (cnet_df[col] != 0).any():
            return True
    return False


def _lonlat_from_campt(
    cnet_df: pd.DataFrame,
    cube_paths: list,
) -> pd.DataFrame:
    """Convert sample/line to lon/lat using ISIS campt via kalasiris.

    For each serial number matched to a cube, writes a coordinate list
    and calls ``campt`` in batch mode to get precise lon/lat from the
    camera model.

    Returns a copy of cnet_df with ``campt_lon`` and ``campt_lat`` columns.

    Raises
    ------
    RuntimeError
        If campt fails for all cubes (e.g. missing ISIS installation).
    """
    import csv
    import os
    import tempfile
    from pathlib import Path

    import kalasiris as isis

    from isistools.io.cubes import build_serial_lookup

    # Use full user environment so ISIS binaries find all required libraries
    isis.environ = os.environ.copy()

    clock_lookup = build_serial_lookup([Path(p) for p in cube_paths])

    df = cnet_df.copy()
    df["campt_lon"] = np.nan
    df["campt_lat"] = np.nan

    n_success = 0
    for sn in df["serialnumber"].unique():
        clock = sn.rsplit("/", 1)[-1]
        if clock not in clock_lookup:
            continue
        cube_path = clock_lookup[clock]
        mask = df["serialnumber"] == sn
        measures = df.loc[mask, ["sample", "line"]]
        if measures.empty:
            continue

        try:
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".csv", delete=False
            ) as coord_f:
                coord_path = coord_f.name
                for _, row in measures.iterrows():
                    coord_f.write(f"{row['sample']},{row['line']}\n")

            out_path = coord_path + ".out.csv"
            isis.campt(
                from_=str(cube_path),
                usecoordlist="true",
                coordlist=coord_path,
                coordtype="image",
                format="flat",
                append="false",
                to=out_path,
            )

            # Parse flat CSV output
            with open(out_path) as f:
                reader = csv.DictReader(f)
                lons = []
                lats = []
                for row in reader:
                    lats.append(float(row["PlanetocentricLatitude"]))
                    lons.append(float(row["PositiveEast360Longitude"]))

            if len(lons) == mask.sum():
                df.loc[mask, "campt_lon"] = lons
                df.loc[mask, "campt_lat"] = lats
                n_success += 1

        except Exception:
            continue
        finally:
            for p in [coord_path, out_path]:
                try:
                    os.unlink(p)
                except OSError:
                    pass

    if n_success == 0:
        raise RuntimeError("campt failed for all cubes")

    return df


def cnet_to_geodataframe(
    cnet_df: pd.DataFrame,
    cube_paths: list | None = None,
) -> gpd.GeoDataFrame:
    """Convert a control network DataFrame to a GeoDataFrame.

    Groups measures by point and creates point geometries. Uses adjusted
    or apriori body-fixed coordinates when available. For pre-jigsaw
    networks without ground coordinates, uses ISIS ``campt`` via
    kalasiris for precise camera-model conversion.

    Parameters
    ----------
    cnet_df : pd.DataFrame
        Control network with one row per measure.
    cube_paths : list, optional
        Cube file paths. Required for footprint-based approximation
        when the network has no ground coordinates.

    Returns
    -------
    gpd.GeoDataFrame
        One row per control point with lon/lat geometry.
    """
    needs_conversion = not _has_ground_coords(cnet_df) and cube_paths is not None

    if needs_conversion:
        cnet_df = _lonlat_from_campt(cnet_df, cube_paths)
        lon_col, lat_col = "campt_lon", "campt_lat"
    else:
        lon_col = None
        lat_col = None
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
        lon = group[lon_col].mean()
        lat = group[lat_col].mean()

        if np.isnan(lon) or np.isnan(lat):
            continue

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

    if not records:
        return gpd.GeoDataFrame(
            columns=["pointId", "geometry", "status", "n_measures",
                     "residual_magnitude", "pointType"],
            geometry="geometry", crs="EPSG:4326",
        )

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
    import holoviews as hv

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
