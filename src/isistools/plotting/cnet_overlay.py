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

    # Rename sample/line to x/y to match the rioxarray image dimensions,
    # so HoloViews shares the same axes when overlaying on an image.
    df = df.rename(columns={"sample": "x", "line": "y"})

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
            x="x",
            y="y",
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
        "registered": {"color": "black", "marker": "circle", "size": 30, "alpha": 0.7},
        "unregistered": {"color": "#e74c3c", "marker": "circle", "size": 40, "alpha": 0.7},
        "ignored": {"color": "red", "marker": "circle", "size": 50, "alpha": 0.9},
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


def _has_lonlat_coords(cnet_df: pd.DataFrame) -> bool:
    """Check if the control network has non-zero lon/lat columns (degrees)."""
    for col in ["adjustedLon", "adjustedLat", "aprioriLon", "aprioriLat"]:
        if col in cnet_df.columns and (cnet_df[col] != 0).any():
            return True
    return False


def _has_bodyfixed_coords(cnet_df: pd.DataFrame) -> bool:
    """Check if the control network has non-zero body-fixed XYZ (meters)."""
    for prefix in ["adjusted", "apriori"]:
        cols = [f"{prefix}X", f"{prefix}Y", f"{prefix}Z"]
        if all(c in cnet_df.columns for c in cols):
            if any((cnet_df[c] != 0).any() for c in cols):
                return True
    return False


def _bodyfixed_to_lonlat(cnet_df: pd.DataFrame) -> tuple[str, str]:
    """Convert body-fixed XYZ columns to lon/lat (degrees).

    ISIS body-fixed coordinates are planetocentric (X, Y, Z) in meters.
    Conversion: lon = atan2(Y, X), lat = atan2(Z, sqrt(X² + Y²)).
    Returns positive-east 0–360 longitude.
    """
    for prefix in ["adjusted", "apriori"]:
        x_col, y_col, z_col = f"{prefix}X", f"{prefix}Y", f"{prefix}Z"
        if all(c in cnet_df.columns for c in [x_col, y_col, z_col]):
            if (cnet_df[x_col] != 0).any():
                lon_col = f"{prefix}Lon"
                lat_col = f"{prefix}Lat"
                cnet_df[lon_col] = np.degrees(np.arctan2(cnet_df[y_col], cnet_df[x_col])) % 360
                cnet_df[lat_col] = np.degrees(
                    np.arctan2(cnet_df[z_col], np.sqrt(cnet_df[x_col]**2 + cnet_df[y_col]**2))
                )
                return lon_col, lat_col
    raise ValueError("No body-fixed XYZ coordinates found")


def _campt_one_serial(cube_path, samples, lines):
    """Run campt for one cube and return (lons, lats) or None on failure."""
    import csv
    import hashlib
    import os
    import tempfile

    import kalasiris as isis

    from isistools.io.cache import get_cache

    # Check cache first — key on cube mtime + coordinate hash
    cache = get_cache()
    coord_bytes = str((samples, lines)).encode()
    coord_hash = hashlib.md5(coord_bytes).hexdigest()[:12]
    try:
        mtime_ns = os.stat(cube_path).st_mtime_ns
    except OSError:
        mtime_ns = 0
    cache_key = f"campt:{cube_path}:{mtime_ns}:{coord_hash}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    isis.environ = os.environ.copy()

    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".csv", delete=False
        ) as coord_f:
            coord_path = coord_f.name
            for s, l in zip(samples, lines):
                coord_f.write(f"{s},{l}\n")

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

        with open(out_path) as f:
            reader = csv.DictReader(f)
            lons = []
            lats = []
            for row in reader:
                lats.append(float(row["PlanetocentricLatitude"]))
                lons.append(float(row["PositiveEast360Longitude"]))

        if len(lons) == len(samples):
            result = (lons, lats)
            cache.set(cache_key, result)
            return result
    except Exception:
        pass
    finally:
        for p in [coord_path, out_path]:
            try:
                os.unlink(p)
            except OSError:
                pass
    return None


def _lonlat_from_campt(
    cnet_df: pd.DataFrame,
    cube_paths: list,
    clock_lookup: dict | None = None,
) -> pd.DataFrame:
    """Convert sample/line to lon/lat using ISIS campt via kalasiris.

    Runs campt calls in parallel using a thread pool (I/O-bound subprocess
    work benefits from threads, not processes).

    Returns a copy of cnet_df with ``campt_lon`` and ``campt_lat`` columns.

    Raises
    ------
    RuntimeError
        If campt fails for all cubes (e.g. missing ISIS installation).
    """
    import os
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from pathlib import Path

    if clock_lookup is None:
        from isistools.io.cubes import build_serial_lookup
        clock_lookup = build_serial_lookup([Path(p) for p in cube_paths])

    df = cnet_df.copy()
    df["campt_lon"] = np.nan
    df["campt_lat"] = np.nan

    # Build work items: (serial_number, cube_path, samples, lines, mask)
    work = []
    for sn in df["serialnumber"].unique():
        clock = sn.rsplit("/", 1)[-1]
        if clock not in clock_lookup:
            continue
        cube_path = clock_lookup[clock]
        mask = df["serialnumber"] == sn
        measures = df.loc[mask, ["sample", "line"]]
        if measures.empty:
            continue
        work.append((sn, cube_path, measures["sample"].tolist(),
                      measures["line"].tolist(), mask))

    n_workers = min(len(work), os.cpu_count() or 4)
    n_success = 0

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(_campt_one_serial, cube_path, samples, lines): (sn, mask)
            for sn, cube_path, samples, lines, mask in work
        }
        for future in as_completed(futures):
            sn, mask = futures[future]
            result = future.result()
            if result is not None:
                lons, lats = result
                df.loc[mask, "campt_lon"] = lons
                df.loc[mask, "campt_lat"] = lats
                n_success += 1

    if n_success == 0:
        raise RuntimeError("campt failed for all cubes")

    return df


def cnet_to_geodataframe(
    cnet_df: pd.DataFrame,
    cube_paths: list | None = None,
    clock_lookup: dict | None = None,
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
        Cube file paths. Required when the network has no ground
        coordinates (pre-jigsaw).
    clock_lookup : dict, optional
        Pre-built clock-count → cube-path mapping. When provided,
        avoids re-parsing cube labels to build the serial lookup.

    Returns
    -------
    gpd.GeoDataFrame
        One row per control point with lon/lat geometry.
    """
    if _has_lonlat_coords(cnet_df):
        # Already have lon/lat in degrees
        lon_col = lat_col = None
        for lon_candidate in ["adjustedLon", "aprioriLon"]:
            if lon_candidate in cnet_df.columns and (cnet_df[lon_candidate] != 0).any():
                lon_col = lon_candidate
                break
        for lat_candidate in ["adjustedLat", "aprioriLat"]:
            if lat_candidate in cnet_df.columns and (cnet_df[lat_candidate] != 0).any():
                lat_col = lat_candidate
                break
    elif _has_bodyfixed_coords(cnet_df):
        # Convert body-fixed XYZ (meters) to lon/lat (degrees)
        cnet_df = cnet_df.copy()
        lon_col, lat_col = _bodyfixed_to_lonlat(cnet_df)
    elif cube_paths is not None:
        # No ground coords at all — use campt to convert sample/line
        cnet_df = _lonlat_from_campt(cnet_df, cube_paths, clock_lookup=clock_lookup)
        lon_col, lat_col = "campt_lon", "campt_lat"
    else:
        raise ValueError(
            "Cannot determine lon/lat for control network points. "
            "No lon/lat columns, no body-fixed XYZ, and no cube_paths for campt. "
            f"Available columns: {list(cnet_df.columns)}"
        )

    # Aggregate per point
    point_groups = cnet_df.groupby("pointId")

    records = []
    for point_id, group in point_groups:
        lon = group[lon_col].mean()
        lat = group[lat_col].mean()

        if np.isnan(lon) or np.isnan(lat) or (lon == 0 and lat == 0):
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

    # Rename to x/y to match rioxarray image dimensions
    df = df.rename(columns={"sample": "x", "line": "y"})

    # Create segments from (x, y) to (x + res*scale, y + res*scale)
    df["x_end"] = df["x"] + df.get("residualSample", 0.0) * scale
    df["y_end"] = df["y"] + df.get("residualLine", 0.0) * scale

    segments = hv.Segments(
        df, kdims=["x", "y", "x_end", "y_end"]
    ).opts(
        color="#e74c3c",
        line_width=1.5,
        alpha=0.8,
    )

    return segments
