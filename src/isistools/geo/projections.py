"""Projection and CRS utilities for ISIS cubes.

Extracts projection information from ISIS cube labels and converts
to pyproj/proj4 representations for use with geopandas and rioxarray.
"""

import pvl
from pyproj import CRS

# Mapping from ISIS projection names to proj4 projection identifiers
ISIS_TO_PROJ4 = {
    "Equirectangular": "eqc",
    "SimpleCylindrical": "eqc",
    "Sinusoidal": "sinu",
    "Mercator": "merc",
    "TransverseMercator": "tmerc",
    "PolarStereographic": "stere",
    "Orthographic": "ortho",
    "LambertConformal": "lcc",
    "LambertAzimuthalEqualArea": "laea",
    "PointPerspective": "nsper",
    "Mollweide": "moll",
    "Robinson": "robin",
}


def mapping_to_crs(mapping: dict) -> CRS:
    """Convert an ISIS Mapping group dict to a pyproj CRS.

    Parameters
    ----------
    mapping : dict
        The Mapping group from an ISIS cube label (as returned by
        ``pvl.load(...)["IsisCube"]["Mapping"]``) or from a MAP file.

    Returns
    -------
    pyproj.CRS
        Coordinate reference system for the projection.
    """
    proj_name = mapping.get("ProjectionName", "Equirectangular")
    proj4_id = ISIS_TO_PROJ4.get(proj_name)
    if proj4_id is None:
        msg = f"Projection '{proj_name}' not supported; add to ISIS_TO_PROJ4"
        raise NotImplementedError(msg)

    # Target body radius
    eq_radius = _to_meters(mapping.get("EquatorialRadius", 3396190.0))
    pol_radius = _to_meters(mapping.get("PolarRadius", eq_radius))

    # Center coordinates
    center_lon = float(mapping.get("CenterLongitude", 0.0))
    center_lat = float(mapping.get("CenterLatitude", 0.0))

    parts = [
        f"+proj={proj4_id}",
        f"+lon_0={center_lon}",
        f"+lat_ts={center_lat}" if proj4_id == "eqc" else f"+lat_0={center_lat}",
        f"+a={eq_radius}",
        f"+b={pol_radius}",
        "+units=m",
        "+no_defs",
        "+type=crs",
    ]

    return CRS.from_proj4(" ".join(parts))


def mapping_to_proj4(mapping: dict) -> str:
    """Convert an ISIS Mapping group dict to a proj4 string.

    Parameters
    ----------
    mapping : dict
        The Mapping group from an ISIS cube label (as returned by
        ``pvl.load(...)["IsisCube"]["Mapping"]``).

    Returns
    -------
    str
        Proj4 projection string.
    """
    return mapping_to_crs(mapping).to_proj4()


def _to_meters(value) -> float:
    """Convert a PVL quantity to meters."""
    if isinstance(value, pvl.Units):
        v = float(value.value)
        unit = str(value.units).lower()
        if unit in ("km", "kilometers"):
            return v * 1000.0
        elif unit in ("m", "meters"):
            return v
        return v  # assume meters
    return float(value)
