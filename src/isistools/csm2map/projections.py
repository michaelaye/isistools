"""Projection and CRS utilities for ISIS cubes.

Extracts projection information from ISIS cube labels and converts
to pyproj/proj4 representations for use with geopandas and rioxarray.

Also provides pure-geometry helpers for converting between the three
ISIS latitude/longitude conventions and the single convention csm2map
uses internally (planetocentric, positive east, 360° domain). ISIS MAP
files may specify any of:

- ``LatitudeType = Planetocentric | Planetographic``
- ``LongitudeDirection = PositiveEast | PositiveWest``
- ``LongitudeDomain = 360 | 180``

and csm2map must normalize all of them before feeding CSM's
``groundToImage`` (which is planetocentric / positive-east) or pyproj's
forward projection (which uses planetocentric when ``+a == +b`` or the
body is spherical). The ``planetographic_to_planetocentric`` and
``planetocentric_to_planetographic`` helpers handle the ~0.3° latitude
shift on Mars; the longitude helpers handle the cosmetic but critical
positive-east/west flip and 180/360 domain rotation.
"""

import math

import numpy as np
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

    Raises
    ------
    ValueError
        If the Mapping group does not carry an ``EquatorialRadius``
        keyword. Previous versions silently defaulted to Mars radii,
        which silently mis-projected any non-Mars MAP file. Callers must
        now either provide explicit body radii in the Mapping group or
        build the CRS themselves with a known body.
    NotImplementedError
        If the projection name is not in :data:`ISIS_TO_PROJ4`.
    """
    proj_name = mapping.get("ProjectionName", "Equirectangular")
    proj4_id = ISIS_TO_PROJ4.get(proj_name)
    if proj4_id is None:
        msg = f"Projection '{proj_name}' not supported; add to ISIS_TO_PROJ4"
        raise NotImplementedError(msg)

    # Target body radius — required; no silent Mars default.
    if "EquatorialRadius" not in mapping:
        target = mapping.get("TargetName", "<unknown>")
        msg = (
            f"Mapping group for target {target!r} is missing EquatorialRadius. "
            f"Cannot build a CRS without explicit body radii. Previous versions "
            f"silently defaulted to Mars; this version refuses to guess. Add "
            f"EquatorialRadius (and PolarRadius) to the MAP file or Mapping "
            f"group, or construct the CRS via TargetBody in csm2map."
        )
        raise ValueError(msg)
    eq_radius = _to_meters(mapping["EquatorialRadius"])
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


def mapping_to_wkt2(mapping: dict) -> str:
    """Convert an ISIS Mapping group dict to a WKT2 CRS string.

    This is the lossless alternative to proj4 strings.

    Parameters
    ----------
    mapping : dict
        The Mapping group from an ISIS cube label (as returned by
        ``pvl.load(...)["IsisCube"]["Mapping"]``).

    Returns
    -------
    str
        WKT2 CRS string.
    """
    return mapping_to_crs(mapping).to_wkt()


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


# ----------------------------------------------------------------------
# Latitude / longitude convention conversions
#
# ISIS MAP files may specify latitude as Planetocentric or Planetographic
# (sometimes spelled "Planetodetic") and longitude as PositiveEast or
# PositiveWest, with a domain of 180° or 360°. csm2map internally works
# in Planetocentric / PositiveEast / 360° because that matches both
# CSM's ECEF geometry and pyproj's spherical projections. These helpers
# convert MAP-file values into csm2map's internal convention.
#
# The planetographic↔planetocentric conversion is a pure geometric
# transform on an oblate spheroid, derived from:
#
#   tan(planetocentric) = (b / a)^2 * tan(planetographic)
#
# where ``a`` is the equatorial (semi-major) radius and ``b`` is the
# polar (semi-minor) radius. The two latitudes agree at the equator
# and at the poles; the maximum difference on Mars (a=3396.19 km,
# b=3376.20 km) is ~0.33° around ±45°. Ignoring the conversion for a
# Mars MAP file with LatitudeType=Planetographic produces a ~20 km
# ground-location error at mid-latitudes — which is why csm2map v0.7.x
# had a silent-bug landmine in this area (the value was read but
# never converted).


def planetographic_to_planetocentric(
    lat_deg: float | np.ndarray,
    eq_radius_m: float,
    polar_radius_m: float,
) -> float | np.ndarray:
    """Convert planetographic (geodetic) latitude to planetocentric.

    Parameters
    ----------
    lat_deg : float or ndarray
        Planetographic latitude(s) in degrees.
    eq_radius_m : float
        Equatorial (semi-major) radius in meters.
    polar_radius_m : float
        Polar (semi-minor) radius in meters.

    Returns
    -------
    float or ndarray
        Planetocentric latitude(s) in degrees, same shape as input.

    Notes
    -----
    For a spherical body (``eq_radius == polar_radius``) this is the
    identity function. Handles the ±90° poles exactly via short-circuit.
    """
    if eq_radius_m == polar_radius_m:
        return lat_deg  # sphere — both conventions are identical

    ratio_sq = (polar_radius_m / eq_radius_m) ** 2

    if np.isscalar(lat_deg):
        if abs(lat_deg) >= 90.0 - 1e-12:
            return float(lat_deg)
        return math.degrees(math.atan(ratio_sq * math.tan(math.radians(lat_deg))))

    lat = np.asarray(lat_deg, dtype=np.float64)
    out = np.empty_like(lat)
    pole_mask = np.abs(lat) >= 90.0 - 1e-12
    out[pole_mask] = lat[pole_mask]
    safe = ~pole_mask
    out[safe] = np.degrees(np.arctan(ratio_sq * np.tan(np.radians(lat[safe]))))
    return out


def planetocentric_to_planetographic(
    lat_deg: float | np.ndarray,
    eq_radius_m: float,
    polar_radius_m: float,
) -> float | np.ndarray:
    """Inverse of :func:`planetographic_to_planetocentric`.

    Uses ``tan(planetographic) = (a / b)^2 * tan(planetocentric)``.
    """
    if eq_radius_m == polar_radius_m:
        return lat_deg

    inv_ratio_sq = (eq_radius_m / polar_radius_m) ** 2

    if np.isscalar(lat_deg):
        if abs(lat_deg) >= 90.0 - 1e-12:
            return float(lat_deg)
        return math.degrees(math.atan(inv_ratio_sq * math.tan(math.radians(lat_deg))))

    lat = np.asarray(lat_deg, dtype=np.float64)
    out = np.empty_like(lat)
    pole_mask = np.abs(lat) >= 90.0 - 1e-12
    out[pole_mask] = lat[pole_mask]
    safe = ~pole_mask
    out[safe] = np.degrees(np.arctan(inv_ratio_sq * np.tan(np.radians(lat[safe]))))
    return out


def normalize_longitude(
    lon_deg: float | np.ndarray,
    *,
    direction: str = "PositiveEast",
    domain: int | str = 360,
) -> float | np.ndarray:
    """Convert a MAP-file longitude value to csm2map's internal convention.

    csm2map internally uses **PositiveEast** longitude in the **360°**
    domain (i.e. lon ∈ [0, 360)). This helper converts any ISIS-legal
    combination to that canonical form.

    Parameters
    ----------
    lon_deg : float or ndarray
        Longitude value(s) in degrees.
    direction : {"PositiveEast", "PositiveWest"}
        The longitude direction of the input values. Case-insensitive.
    domain : {180, 360, "180", "360"}
        The longitude domain of the input values. ``360`` means the
        input is in ``[0, 360)``; ``180`` means ``[-180, 180)``.

    Returns
    -------
    float or ndarray
        Longitude values in PositiveEast, 360° domain.

    Notes
    -----
    - Positive-west → positive-east: ``lon_pe = -lon_pw`` (mod 360).
    - 180° domain → 360° domain: ``lon_360 = lon_180 % 360`` (handles
      negative values by wrapping).
    - Output is always in ``[0, 360)`` regardless of input range.
    """
    direction_norm = str(direction).strip().lower().replace(" ", "")
    if direction_norm in ("positivewest", "pw", "west"):
        lon = -np.asarray(lon_deg, dtype=np.float64)
    elif direction_norm in ("positiveeast", "pe", "east"):
        lon = np.asarray(lon_deg, dtype=np.float64)
    else:
        msg = f"Unrecognized longitude direction {direction!r}"
        raise ValueError(msg)

    # Normalize to [0, 360) domain regardless of the stated input domain
    # (this also handles negative values from a 180° domain or from the
    # positive-west flip above).
    lon = np.mod(lon, 360.0)

    if np.isscalar(lon_deg):
        return float(lon)
    return lon


def normalize_latitude_from_mapping(
    lat_deg: float | np.ndarray,
    mapping: dict,
    eq_radius_m: float,
    polar_radius_m: float,
) -> float | np.ndarray:
    """Convert a latitude read from a Mapping group to planetocentric.

    Reads the ``LatitudeType`` keyword from the Mapping group. If absent
    or ``Planetocentric``, the input is returned unchanged. If
    ``Planetographic`` (or the synonym ``Planetodetic``), the value is
    converted via :func:`planetographic_to_planetocentric`.

    Any other string raises :class:`ValueError`.
    """
    raw = mapping.get("LatitudeType", "Planetocentric")
    lat_type = str(raw).strip().lower()
    if lat_type in ("planetocentric", "centric"):
        return lat_deg
    if lat_type in ("planetographic", "planetodetic", "geodetic", "graphic"):
        return planetographic_to_planetocentric(lat_deg, eq_radius_m, polar_radius_m)
    msg = f"Unrecognized LatitudeType {raw!r}; expected Planetocentric or Planetographic"
    raise ValueError(msg)


def normalize_longitude_from_mapping(
    lon_deg: float | np.ndarray,
    mapping: dict,
) -> float | np.ndarray:
    """Convert a longitude read from a Mapping group to PositiveEast / 360°.

    Reads ``LongitudeDirection`` (default ``PositiveEast``) and
    ``LongitudeDomain`` (default ``360``) from the Mapping group and
    calls :func:`normalize_longitude`.
    """
    direction = mapping.get("LongitudeDirection", "PositiveEast")
    domain = mapping.get("LongitudeDomain", 360)
    return normalize_longitude(lon_deg, direction=direction, domain=domain)
