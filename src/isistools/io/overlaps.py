"""Parse findimageoverlaps output into GeoDataFrames.

The ``overlap_list.lis`` produced by ISIS ``findimageoverlaps`` alternates:

- **Odd lines:** comma-separated serial numbers (the images involved)
- **Even lines:** WKB hex encoding of the overlap polygon

Entries with a single serial number are individual image footprints.
Entries with 2+ serial numbers are overlap zones.
"""

from __future__ import annotations

import binascii
from pathlib import Path

import geopandas as gpd
from shapely import wkb


def parse_overlap_list(path: str | Path) -> gpd.GeoDataFrame:
    """Parse an ISIS overlap_list.lis into a GeoDataFrame.

    Parameters
    ----------
    path : path-like
        Path to the overlap list file produced by ``findimageoverlaps``.

    Returns
    -------
    gpd.GeoDataFrame
        One row per overlap zone with columns: ``serials`` (raw string),
        ``images`` (list of serial numbers), ``n_images`` (count),
        ``zone_type`` (footprint / 2-way overlap / N-way overlap),
        ``area_deg2`` (polygon area in square degrees), and ``geometry``.
    """
    path = Path(path)
    with open(path) as f:
        lines = [line.strip() for line in f if line.strip()]

    records = []
    for i in range(0, len(lines), 2):
        serials_raw = lines[i]
        hex_data = lines[i + 1]

        serial_list = [s.strip() for s in serials_raw.split(",")]
        geom = wkb.loads(binascii.unhexlify(hex_data))

        n_images = len(serial_list)
        if n_images == 1:
            zone_type = "footprint"
        elif n_images == 2:
            zone_type = "2-way overlap"
        else:
            zone_type = f"{n_images}-way overlap"

        records.append(
            {
                "serials": serials_raw,
                "images": serial_list,
                "n_images": n_images,
                "zone_type": zone_type,
                "area_deg2": geom.area,
                "geometry": geom,
            }
        )

    gdf = gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")
    return gdf
