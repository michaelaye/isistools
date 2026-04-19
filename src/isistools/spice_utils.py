"""SPICE utilities for isistools — Sun distance and ephemeris queries.

Uses spiceypy for direct SPICE kernel queries. This is the right tool for
simple ephemeris lookups (Sun distance, body positions). For camera model
construction, use ALE (ale.loads) instead.

Requires spiceypy and generic SPICE kernels from ISISDATA:
  - base/kernels/lsk/naif????.tls  (leap second)
  - base/kernels/spk/de???.bsp    (planetary ephemeris)
  - base/kernels/spk/mar???.bsp   (Mars satellites)
  - base/kernels/pck/pck?????.tpc (planetary constants)
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np

_ISISDATA = Path(os.environ.get("ISISDATA", str(Path.home() / "Dropbox" / "data" / "isisdata")))


def _require_spiceypy():
    try:
        import spiceypy

        return spiceypy
    except ImportError:
        raise ImportError(
            "spiceypy is required for SPICE computations. "
            "Install with: conda install -c conda-forge spiceypy"
        )


def sun_distance_au(
    utc_time: str,
    target: str = "mars",
    isisdata: Path | None = None,
) -> float:
    """Compute Sun-to-target distance in AU at a given UTC time.

    Parameters
    ----------
    utc_time : str
        UTC time string, e.g. "2008-12-21T10:01:34.159".
    target : str
        Target body name (default "mars").
    isisdata : Path, optional
        ISISDATA directory. Defaults to $ISISDATA or ~/Dropbox/data/isisdata.

    Returns
    -------
    float
        Distance in AU.
    """
    dist_km = sun_distance_km(utc_time, target=target, isisdata=isisdata)
    return dist_km / 1.49597870691e8


def sun_distance_km(
    utc_time: str,
    target: str = "mars",
    isisdata: Path | None = None,
) -> float:
    """Compute Sun-to-target distance in km at a given UTC time.

    Parameters
    ----------
    utc_time : str
        UTC time string, e.g. "2008-12-21T10:01:34.159".
    target : str
        Target body name (default "mars").
    isisdata : Path, optional
        ISISDATA directory.

    Returns
    -------
    float
        Distance in km.
    """
    spice = _require_spiceypy()
    if isisdata is None:
        isisdata = _ISISDATA
    base_kernels = isisdata / "base" / "kernels"

    lsk = sorted(base_kernels.glob("lsk/naif????.tls"))[-1]
    spk_de = sorted(base_kernels.glob("spk/de???.bsp"))[-1]
    spk_mar = sorted(base_kernels.glob("spk/mar???.bsp"))[-1]
    pck = sorted(base_kernels.glob("pck/pck?????.tpc"))[-1]

    loaded = []
    try:
        for k in [lsk, spk_de, spk_mar, pck]:
            spice.furnsh(str(k))
            loaded.append(str(k))

        # Clean UTC string (remove timezone suffixes pvl may add)
        utc_clean = utc_time.replace("+00:00", "").replace("Z", "")
        et = spice.utc2et(utc_clean)

        # ISIS uses spkpos(target, et, "J2000", "LT+S", "sun") for HiRISE
        # and spkezr("sun", et, "iau_mars", "LT+S", "mars") for CTX
        # Both give the same distance (just different frames/directions)
        sunpos, lt = spice.spkpos(target, et, "J2000", "LT+S", "sun")
        dist = np.linalg.norm(sunpos[:3])
    finally:
        for k in loaded:
            spice.unload(k)

    return float(dist)


def sun_distance_from_cube(cube_path: str | Path) -> float:
    """Compute Sun-to-target distance in km from an ISIS cube label.

    Reads StartTime and TargetName from the cube label.

    Parameters
    ----------
    cube_path : path-like
        Path to a spiceinit'd ISIS cube.

    Returns
    -------
    float
        Distance in km.
    """
    import pvl

    label = pvl.load(str(cube_path))
    inst = label["IsisCube"]["Instrument"]

    start_time_raw = inst["StartTime"]
    if hasattr(start_time_raw, "strftime"):
        start_time = start_time_raw.strftime("%Y-%m-%dT%H:%M:%S.%f")
    else:
        start_time = str(start_time_raw)

    target = str(inst.get("TargetName", "Mars"))
    # Normalize target names
    target_lower = target.lower()
    if target_lower in ("sky", "cal", "phobos", "deimos"):
        target = "Mars"

    return sun_distance_km(start_time, target=target)
