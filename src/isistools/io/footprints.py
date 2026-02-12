"""Read footprint polygons from ISIS cube files.

After running ISIS's ``footprintinit``, each cube contains a polygon
blob that encodes the image footprint on the target body surface.
This module reads those blobs and returns them as Shapely geometries
collected into a GeoDataFrame.

The polygon is stored as a serialized blob in the cube file. The label
contains a pointer (``^Polygon``) giving the byte offset and optionally
the format. The blob is typically GML or WKT text.
"""

from __future__ import annotations

import struct
from pathlib import Path
from typing import TYPE_CHECKING

import geopandas as gpd
import numpy as np
import pandas as pd
import pvl
import shapely
from shapely import wkt as shapely_wkt

if TYPE_CHECKING:
    pass


class FootprintNotFoundError(Exception):
    """Raised when a cube has no footprint polygon.

    This usually means ``footprintinit`` has not been run on the cube.
    """


def _find_polygon_blob(cube_path: Path, label: pvl.PVLModule) -> str:
    """Read the raw polygon blob text from an ISIS cube file.

    ISIS stores the polygon as a text blob (GML or WKT) at a byte
    offset indicated by the ``^Polygon`` pointer in the label.
    The blob is preceded by a small header with the object name
    and size.

    Parameters
    ----------
    cube_path : Path
        Path to the ISIS cube file.
    label : pvl.PVLModule
        Pre-parsed PVL label.

    Returns
    -------
    str
        Raw polygon text (GML or WKT).
    """
    # ISIS stores the footprint as an ``Object = Polygon`` in the label
    # with ``StartByte`` and ``Bytes`` fields indicating where the blob
    # lives in the file.
    polygon_obj = label.get("Polygon")
    if polygon_obj is None:
        raise FootprintNotFoundError(
            f"No Polygon object found in {cube_path}. "
            "Run footprintinit first."
        )

    start_byte = int(polygon_obj["StartByte"]) - 1  # PVL is 1-based
    nbytes = int(polygon_obj["Bytes"])

    with open(cube_path, "rb") as f:
        f.seek(start_byte)
        chunk = f.read(nbytes)

    return chunk.decode("ascii", errors="replace").strip().rstrip("\x00")


def _extract_wkt(text: str) -> str:
    """Extract a complete WKT geometry string from text.

    Handles balanced parentheses to find the end of the WKT.
    """
    depth = 0
    end = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    return text[:end] if end > 0 else text.split("\x00")[0].strip()


def _parse_polygon_text(text: str) -> shapely.Geometry:
    """Parse polygon text (WKT or GML) into a Shapely geometry."""
    text = text.strip()

    # Try WKT first
    if text.startswith(("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION")):
        return shapely_wkt.loads(text)

    # Try GML
    try:
        from shapely import gml
        return gml.loads(text)
    except (ImportError, Exception):
        pass

    # Try ogr as fallback for GML
    try:
        from osgeo import ogr
        geom = ogr.CreateGeometryFromGML(text)
        if geom is not None:
            return shapely_wkt.loads(geom.ExportToWkt())
    except ImportError:
        pass

    raise FootprintNotFoundError(
        f"Could not parse polygon text. First 200 chars: {text[:200]}"
    )


def read_footprint(
    cube_path: str | Path,
    label: pvl.PVLModule | None = None,
) -> shapely.Geometry:
    """Read the footprint polygon from an ISIS cube.

    The cube must have had ``footprintinit`` run on it.

    Parameters
    ----------
    cube_path : path-like
        Path to the ISIS cube file.
    label : pvl.PVLModule, optional
        Pre-parsed PVL label. Parsed from *cube_path* if not given.

    Returns
    -------
    shapely.Geometry
        Footprint polygon in target body lon/lat coordinates.

    Raises
    ------
    FootprintNotFoundError
        If no polygon blob is found (footprintinit not run).
    """
    cube_path = Path(cube_path)
    if label is None:
        label = pvl.load(str(cube_path))
    blob_text = _find_polygon_blob(cube_path, label)
    return _parse_polygon_text(blob_text)


def read_cube_list(cube_list_path: str | Path) -> list[Path]:
    """Read an ISIS-style cube list file (one path per line).

    Skips blank lines and lines starting with ``#``.
    """
    cube_list_path = Path(cube_list_path)
    paths = []
    with open(cube_list_path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                paths.append(Path(line))
    return paths


def load_footprints(
    cube_list: str | Path | list[str | Path],
    skip_errors: bool = False,
) -> gpd.GeoDataFrame:
    """Load footprints from all cubes in a list.

    Parameters
    ----------
    cube_list : path-like or list
        Either a path to an ISIS cube list file, or a list of
        individual cube file paths.
    skip_errors : bool
        If True, cubes without footprints are silently skipped.
        If False, raises on the first failure.

    Returns
    -------
    gpd.GeoDataFrame
        Footprints with columns: path, filename, geometry,
        and selected label metadata.
    """
    if isinstance(cube_list, (str, Path)):
        cube_list_path = Path(cube_list)
        if cube_list_path.is_file() and cube_list_path.suffix in (".lis", ".txt", ".list", ""):
            cube_paths = read_cube_list(cube_list_path)
        else:
            # Single cube file
            cube_paths = [cube_list_path]
    else:
        cube_paths = [Path(p) for p in cube_list]

    records = []
    for cp in cube_paths:
        try:
            label = pvl.load(str(cp))
            geom = read_footprint(cp, label=label)

            # Extract metadata
            inst = label["IsisCube"].get("Instrument", {})
            mapping = label["IsisCube"].get("Mapping", {})

            clock = str(
                inst.get("SpacecraftClockCount",
                         inst.get("SpacecraftClockStartCount", ""))
            )

            records.append({
                "path": str(cp),
                "filename": cp.name,
                "geometry": geom,
                "target": mapping.get(
                    "TargetName",
                    inst.get("TargetName", "Unknown"),
                ),
                "start_time": str(inst.get("StartTime", "")),
                "instrument": inst.get("InstrumentId", "Unknown"),
                "spacecraft": inst.get(
                    "SpacecraftName",
                    inst.get("SpacecraftId", "Unknown"),
                ),
                "clock": clock,
                "level": 2 if mapping else 1,
            })
        except FootprintNotFoundError:
            if skip_errors:
                continue
            raise
        except Exception as e:
            if skip_errors:
                continue
            raise RuntimeError(f"Failed to read footprint from {cp}") from e

    if not records:
        return gpd.GeoDataFrame(
            columns=["path", "filename", "geometry", "target",
                      "start_time", "instrument", "spacecraft", "clock", "level"],
            geometry="geometry",
        )

    gdf = gpd.GeoDataFrame(records, geometry="geometry")

    # Set CRS — ISIS footprints are in target body lon/lat (planetocentric).
    # We use EPSG:4326 as a proxy; for non-Earth bodies this is approximate
    # but sufficient for visualization.
    gdf = gdf.set_crs(epsg=4326)

    return gdf
