"""Load CSM camera models from ISIS cubes via ale.

This module wraps ale + usgscsm to produce a CSM sensor model from a
spiceinit'd ISIS cube.  The resulting model supports groundToImage() and
imageToGround() calls needed for map projection.
"""

import ctypes
import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import csmapi

_CSM_PLUGIN_LOADED = False


def _ensure_csm_plugin_loaded() -> None:
    """Load the usgscsm CSM plugin via ctypes.

    The conda-forge ``usgscsm`` package installs the plugin as
    ``lib/csmplugins/libusgscsm.dylib`` (or ``.so`` on Linux) but does not
    provide a Python module. We load it manually so the plugin registers
    itself with csmapi on first use.
    """
    global _CSM_PLUGIN_LOADED
    if _CSM_PLUGIN_LOADED:
        return

    # Try explicit plugin path first (conda-forge layout)
    candidates = []
    env_prefix = os.environ.get("CONDA_PREFIX") or os.environ.get("PREFIX")
    if env_prefix:
        candidates.append(Path(env_prefix) / "lib" / "csmplugins" / "libusgscsm.dylib")
        candidates.append(Path(env_prefix) / "lib" / "csmplugins" / "libusgscsm.so")

    # Also check sys.prefix
    import sys

    candidates.append(Path(sys.prefix) / "lib" / "csmplugins" / "libusgscsm.dylib")
    candidates.append(Path(sys.prefix) / "lib" / "csmplugins" / "libusgscsm.so")

    for lib_path in candidates:
        if lib_path.exists():
            ctypes.CDLL(str(lib_path))
            _CSM_PLUGIN_LOADED = True
            return

    # Fallback: try library search path (may work if CSM_PLUGIN_PATH is set)
    for name in ("libusgscsm.dylib", "libusgscsm.so"):
        try:
            ctypes.CDLL(name)
            _CSM_PLUGIN_LOADED = True
            return
        except OSError:
            continue

    msg = (
        "Could not locate libusgscsm plugin. Expected at "
        "$CONDA_PREFIX/lib/csmplugins/libusgscsm.dylib"
    )
    raise RuntimeError(msg)


def load_camera(
    cube_path: str | Path,
    spice_source: str = "isis",
    refresh_isd: bool = True,
) -> "csmapi.RasterGM":
    """Load a CSM RasterGM sensor model from a spiceinit'd ISIS cube.

    Uses ale to generate an ISD (Instrument Support Data) from the cube,
    then iterates usgscsm plugins to construct a CSM model. This replaces
    the knoten.csm.create_csm() call to avoid the knoten dependency
    (which has a broken csmapi conda dep).

    Parameters
    ----------
    cube_path : path-like
        Path to a Level 1 ISIS cube that has been run through spiceinit.
    spice_source : {"isis", "naif", "auto"}
        Where ALE should source SPICE data from when generating the ISD:

        - ``"isis"`` (default): use the SPICE blobs **embedded in the cube**
          itself (``IsisSpice`` driver path). This is the right choice for
          any pipeline that runs ``jigsaw update=true``, because jigsaw
          updates the cube's embedded blobs but does NOT update the live
          NAIF kernels — so reading the live kernels would silently
          discard the jigsaw bundle adjustment results.
        - ``"naif"``: force the live NAIF kernel path (``NaifSpice``
          driver). Use this only when you know the cube has not been
          jigsaw-updated, or when you want to compare against the
          pristine pre-jigsaw geometry.
        - ``"auto"``: let ALE pick the driver via its default heuristic,
          which currently prefers NAIF over ISIS regardless of whether
          the cube has been jigsaw-updated. **Not recommended** for
          jigsaw pipelines.
    refresh_isd : bool
        If True (default), regenerate the ISD JSON every time. Set to
        False to reuse a cached ISD if one already exists alongside the
        cube. Defaults to True because a stale ISD from a pre-jigsaw run
        is the exact failure mode we want to prevent.

    Returns
    -------
    csmapi.RasterGM
        A CSM raster ground-to-image sensor model.

    Raises
    ------
    RuntimeError
        If ale cannot generate an ISD or usgscsm cannot build the model.
    """
    import ale
    import csmapi

    _ensure_csm_plugin_loaded()

    cube_path = Path(cube_path)
    cube_str = str(cube_path)

    # Generate ISD JSON via ALE.  We default to ``only_isis_spice=True``
    # so that the cube's embedded SPICE blobs (which jigsaw updates with
    # update=true) are the source of truth, NOT the live NAIF kernels.
    # csmapi.Isd() reads the JSON file with the same stem as the cube
    # (e.g. ``foo.cub`` -> ``foo.json``), so ale's output must live
    # alongside the cube.
    isd_path = cube_path.with_suffix(".json")
    if refresh_isd or not isd_path.exists():
        if spice_source == "isis":
            isd_str = ale.loads(cube_str, only_isis_spice=True)
        elif spice_source == "naif":
            isd_str = ale.loads(cube_str, only_naif_spice=True)
        elif spice_source == "auto":
            isd_str = ale.loads(cube_str)
        else:
            msg = f"Invalid spice_source={spice_source!r}; must be 'isis', 'naif', or 'auto'"
            raise ValueError(msg)
        isd_path.write_text(isd_str)

    isd = csmapi.Isd(cube_str)

    # Iterate registered CSM plugins and try to construct a model
    for plugin in csmapi.Plugin.getList():
        for model_index in range(plugin.getNumModels()):
            model_name = plugin.getModelName(model_index)
            warnings = csmapi.WarningList()
            if plugin.canModelBeConstructedFromISD(isd, model_name, warnings):
                model = plugin.constructModelFromISD(isd, model_name)
                if isinstance(model, csmapi.RasterGM):
                    return model

    msg = f"No CSM plugin could construct a model from {cube_path}"
    raise RuntimeError(msg)


def get_image_size(model: "csmapi.RasterGM") -> tuple[int, int]:
    """Return (n_lines, n_samples) from a CSM model."""
    size = model.getImageSize()
    return int(size.line), int(size.samp)


def get_target_radii(cube_path: str | Path) -> tuple[float, float]:
    """Extract target body equatorial and polar radii from cube labels.

    Returns
    -------
    tuple of (equatorial_radius_m, polar_radius_m)
    """
    from isistools.io.cubes import read_label

    label = read_label(cube_path)
    # After spiceinit, radii are in the NaifKeywords group
    try:
        naif = label["NaifKeywords"]
        # BODY_RADII is [a, b, c] in km
        radii = naif["BODY499_RADII"]  # Mars = 499; TODO: generalize
        eq_r = float(radii[0]) * 1000.0
        polar_r = float(radii[2]) * 1000.0
    except (KeyError, IndexError):
        # Fallback: Mars defaults
        eq_r = 3396190.0
        polar_r = 3376200.0

    return eq_r, polar_r


def ground_to_image_batch(
    model: "csmapi.RasterGM",
    lats: np.ndarray,
    lons: np.ndarray,
    radii: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Batch groundToImage: lat/lon/radius arrays -> input (lines, samples).

    Parameters
    ----------
    model : csmapi.RasterGM
        CSM sensor model.
    lats, lons : ndarray
        Planetocentric latitude and longitude in radians.
    radii : ndarray
        Surface radius in meters at each point (sphere or DEM).

    Returns
    -------
    lines, samples : ndarray
        Input image coordinates.  Invalid points get NaN.
    """
    import csmapi

    flat = lats.ravel()
    flon = lons.ravel()
    frad = radii.ravel()

    # Convert planetocentric lat/lon/radius to body-fixed ECEF (X, Y, Z)
    cos_lat = np.cos(flat)
    sin_lat = np.sin(flat)
    cos_lon = np.cos(flon)
    sin_lon = np.sin(flon)

    x = frad * cos_lat * cos_lon
    y = frad * cos_lat * sin_lon
    z = frad * sin_lat

    out_lines = np.full(flat.shape, np.nan)
    out_samps = np.full(flat.shape, np.nan)

    # TODO: This loop is the bottleneck to optimize later.
    # Options: Cython wrapper, ctypes batch call, or multiprocessing.
    for i in range(len(flat)):
        if np.isnan(flat[i]):
            continue
        try:
            gp = csmapi.EcefCoord(x[i], y[i], z[i])
            ip = model.groundToImage(gp)
            out_lines[i] = ip.line
            out_samps[i] = ip.samp
        except Exception:
            pass  # leave as NaN

    return (
        out_lines.reshape(lats.shape),
        out_samps.reshape(lats.shape),
    )
