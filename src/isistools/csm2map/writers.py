"""Output writers for map-projected data."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import rasterio

from isistools.csm2map.grid import OutputGrid

if TYPE_CHECKING:
    from isistools.csm2map.camera import TargetBody


def write_geotiff(
    output_path: str | Path,
    data: np.ndarray,
    grid: OutputGrid,
    nodata: float = 0.0,
) -> Path:
    """Write map-projected data as a GeoTIFF.

    Parameters
    ----------
    output_path : path-like
        Output file path (should end in .tif).
    data : ndarray, shape (height, width) or (n_bands, height, width)
        Map-projected pixel data.
    grid : OutputGrid
        Grid definition with CRS and affine transform.
    nodata : float
        NoData value.

    Returns
    -------
    Path to written file.
    """
    output_path = Path(output_path)

    if data.ndim == 2:
        data = data[np.newaxis, ...]
    n_bands, height, width = data.shape

    assert height == grid.height, f"Data height {height} != grid height {grid.height}"
    assert width == grid.width, f"Data width {width} != grid width {grid.width}"

    # Replace NaN with nodata
    data = np.where(np.isnan(data), nodata, data)

    # ZSTD compression is dramatically faster than LZW at similar or
    # better compression ratios for float32 image data. ZSTD levels
    # 1-3 complete in well under half the wall time of LZW on a
    # 100 MB float32 cube while producing slightly smaller output.
    # The NUM_THREADS=ALL_CPUS option lets GDAL parallelize the
    # compression across cores.
    profile = {
        "driver": "GTiff",
        "dtype": data.dtype,
        "width": width,
        "height": height,
        "count": n_bands,
        "crs": grid.crs.to_wkt(),
        "transform": grid.transform,
        "nodata": nodata,
        "compress": "zstd",
        "zstd_level": 3,
        "num_threads": "ALL_CPUS",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(str(output_path), "w", **profile) as dst:
        for b in range(n_bands):
            dst.write(data[b], b + 1)

    return output_path


# Map from EPSG/WKT2 method names to ISIS projection names.
_METHOD_TO_ISIS = {
    "Equidistant Cylindrical": "Equirectangular",
    "Equidistant Cylindrical (Spherical)": "Equirectangular",
    "Sinusoidal": "Sinusoidal",
    "Mercator (variant A)": "Mercator",
    "Mercator (variant B)": "Mercator",
    "Transverse Mercator": "TransverseMercator",
    "Polar Stereographic (variant A)": "PolarStereographic",
    "Polar Stereographic (variant B)": "PolarStereographic",
    "Orthographic": "Orthographic",
    "Lambert Conic Conformal (2SP)": "LambertConformal",
    "Lambert Azimuthal Equal Area": "LambertAzimuthalEqualArea",
    "Lambert Azimuthal Equal Area (Spherical)": "LambertAzimuthalEqualArea",
}

# Legacy: proj4 identifiers (fallback)
_PROJ4_TO_ISIS = {
    "eqc": "Equirectangular",
    "sinu": "Sinusoidal",
    "merc": "Mercator",
    "tmerc": "TransverseMercator",
    "stere": "PolarStereographic",
    "ortho": "Orthographic",
    "lcc": "LambertConformal",
    "laea": "LambertAzimuthalEqualArea",
    "nsper": "PointPerspective",
    "moll": "Mollweide",
    "robin": "Robinson",
}


def write_mapping_pvl(
    output_path: Path,
    grid: OutputGrid,
    body: "TargetBody",
) -> Path:
    """Write an ISIS-compatible Mapping PVL sidecar next to a GeoTIFF.

    The sidecar carries the same metadata that ISIS ``cam2map`` writes
    into a projected cube's ``IsisCube.Mapping`` group: projection name,
    body radii, lat/lon type and direction, ground range, pixel
    resolution, and the ``UpperLeftCornerX/Y`` that pin the grid origin.
    This makes the GeoTIFF interoperable with ISIS workflows that expect
    a PVL Mapping group — e.g. ``automos``, ``mapmos``, or scripts that
    parse ``catlab`` output.

    Parameters
    ----------
    output_path : Path
        Path to the GeoTIFF whose sidecar this is. The PVL is written to
        ``output_path.with_suffix('.pvl')``.
    grid : OutputGrid
        Grid definition (CRS, affine, dimensions, ground range).
    body : TargetBody
        Target body identity and ellipsoid.

    Returns
    -------
    Path to the written ``.pvl`` file.
    """
    # Extract projection name and parameters from the CRS using pyproj's
    # coordinate_operation API (lossless, unlike to_proj4).
    center_lon = 0.0
    center_lat = 0.0
    op = grid.crs.coordinate_operation
    if op is not None:
        method_name = op.method_name
        isis_proj_name = _METHOD_TO_ISIS.get(method_name, method_name)
        # Extract center lon/lat from operation parameters
        for p in op.params:
            name_lower = p.name.lower()
            if "longitude" in name_lower and "origin" in name_lower:
                center_lon = p.value
            elif "standard parallel" in name_lower:
                center_lat = p.value
            elif "latitude" in name_lower and "origin" in name_lower:
                if center_lat == 0.0:  # don't override lat_ts
                    center_lat = p.value
    else:
        isis_proj_name = "Unknown"

    pvl_path = output_path.with_suffix(".pvl")

    # Normalize longitudes to the declared domain [0, 360) so the PVL
    # is self-consistent. _derive_ground_range may return values outside
    # [-180, +180] for antimeridian-crossing strips (e.g. lon_min=179.5,
    # lon_max=181.5); ISIS tools expect the values to match the domain.
    lon_min_pvl = grid.lon_min % 360.0
    lon_max_pvl = grid.lon_max % 360.0
    center_lon_pvl = center_lon % 360.0
    # If wraparound made max < min (e.g. 359° to 1°), that's correct
    # for the 360 domain — ISIS interprets it as crossing 0°.

    lines = [
        "Group = Mapping",
        f"  ProjectionName     = {isis_proj_name}",
        f"  TargetName         = {body.name}",
        f"  EquatorialRadius   = {body.radius_equatorial_m:.1f} <meters>",
        f"  PolarRadius        = {body.radius_polar_m:.1f} <meters>",
        "  LatitudeType       = Planetocentric",
        "  LongitudeDirection = PositiveEast",
        "  LongitudeDomain    = 360",
        f"  CenterLatitude     = {center_lat}",
        f"  CenterLongitude    = {center_lon_pvl}",
        f"  MinimumLatitude    = {grid.lat_min}",
        f"  MaximumLatitude    = {grid.lat_max}",
        f"  MinimumLongitude   = {lon_min_pvl}",
        f"  MaximumLongitude   = {lon_max_pvl}",
        f"  PixelResolution    = {grid.resolution} <meters/pixel>",
        f"  UpperLeftCornerX   = {grid.transform.c} <meters>",
        f"  UpperLeftCornerY   = {grid.transform.f} <meters>",
        "End_Group",
        "",
        "Group = Dimensions",
        f"  Samples            = {grid.width}",
        f"  Lines              = {grid.height}",
        "End_Group",
        "",
        "Group = AlgorithmName",
        "  Name               = csm2map",
        "  Version            = isistools",
        "End_Group",
    ]

    pvl_path.write_text("\n".join(lines) + "\n")
    return pvl_path
