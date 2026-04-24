"""Load CSM camera models from ISIS cubes via ale.

This module wraps ale + usgscsm to produce a CSM sensor model from a
spiceinit'd ISIS cube.  The resulting model supports groundToImage() and
imageToGround() calls needed for map projection.

Alongside the CSM model, ``load_camera`` also returns a :class:`TargetBody`
describing the target body's ellipsoid (semi-axes, mean radius, NAIF ID).
This lets csm2map project cubes from any mission without hardcoding Mars-
specific constants — the body info comes straight from ALE's ISD, which
is computed from the same SPICE kernels the camera model itself uses.
"""

import ctypes
import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    import csmapi

_CSM_PLUGIN_LOADED = False


@dataclass(frozen=True)
class TargetBody:
    """Target body ellipsoid + identity, parsed from an ALE ISD.

    All radii are stored in **meters** regardless of the ISD unit, so
    downstream code never has to unit-convert.

    Attributes
    ----------
    name : str
        Human-readable target name (e.g. ``"MARS"``, ``"EUROPA"``). Taken
        from the cube label's ``IsisCube.Instrument.TargetName`` when
        available, otherwise synthesized from the NAIF body code as
        ``"BODY<naif_id>"``.
    naif_id : int
        NAIF body code (e.g. 499 = Mars, 301 = Moon, 502 = Europa). Read
        directly from ``isd["naif_keywords"]["BODY_CODE"]``.
    radius_equatorial_m : float
        Equatorial (semi-major) radius in meters.
    radius_polar_m : float
        Polar (semi-minor) radius in meters.
    radius_mean_m : float
        ``(2 * radius_equatorial_m + radius_polar_m) / 3`` — the IAU
        convention for a volume-equivalent spherical radius. Used as the
        fallback radius for DEM lookups where the DEM reports nodata.
    """

    name: str
    naif_id: int
    radius_equatorial_m: float
    radius_polar_m: float
    radius_mean_m: float

    @classmethod
    def from_isd(cls, isd: dict[str, Any], target_name: str | None = None) -> "TargetBody":
        """Construct a :class:`TargetBody` from an ALE ISD dict.

        Parameters
        ----------
        isd : dict
            Parsed ALE ISD (e.g. ``json.loads(ale.loads(cube))``). Must
            contain a top-level ``radii`` dict with ``semimajor``,
            ``semiminor``, ``unit`` keys, plus a ``naif_keywords`` dict
            carrying ``BODY_CODE``.
        target_name : str, optional
            Human-readable target name (normally the cube label's
            ``Instrument.TargetName``). If omitted the body name is
            synthesized from the NAIF ID as ``"BODY<naif_id>"``.

        Raises
        ------
        KeyError
            If required ISD fields are missing.
        ValueError
            If the ISD's ``radii`` dict uses an unrecognized unit, or if
            the cube's ``NaifKeywords.BODY<code>_RADII`` cross-check
            disagrees with ``isd["radii"]`` by more than 1 meter. This
            catches stale SPICE blobs and corrupted cubes.
        """
        radii_block = isd["radii"]
        unit = str(radii_block.get("unit", "km")).lower()
        if unit in ("km", "kilometers"):
            scale = 1000.0
        elif unit in ("m", "meters"):
            scale = 1.0
        else:
            msg = f"Unrecognized radii unit in ISD: {unit!r} (expected 'km' or 'm')"
            raise ValueError(msg)

        eq_m = float(radii_block["semimajor"]) * scale
        polar_m = float(radii_block["semiminor"]) * scale

        naif_kw = isd.get("naif_keywords", {})
        naif_id = int(naif_kw["BODY_CODE"])

        # Cross-check against the explicit BODY<code>_RADII keyword if
        # present. Both values come from the same SPICE kernel load so
        # they should agree to floating-point precision; a mismatch
        # indicates corrupted blobs or a kernel mixup.
        body_radii_key = f"BODY{naif_id}_RADII"
        if body_radii_key in naif_kw:
            br = naif_kw[body_radii_key]
            br_eq_m = float(br[0]) * scale
            br_polar_m = float(br[2]) * scale
            if abs(br_eq_m - eq_m) > 1.0 or abs(br_polar_m - polar_m) > 1.0:
                msg = (
                    f"ISD radii {radii_block} disagree with {body_radii_key}={br} "
                    f"(unit={unit}). Stale SPICE blob or corrupted cube?"
                )
                raise ValueError(msg)

        name = (target_name or f"BODY{naif_id}").upper()
        return cls(
            name=name,
            naif_id=naif_id,
            radius_equatorial_m=eq_m,
            radius_polar_m=polar_m,
            radius_mean_m=(2.0 * eq_m + polar_m) / 3.0,
        )


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
) -> tuple["csmapi.RasterGM", TargetBody]:
    """Load a CSM RasterGM sensor model + target body info from an ISIS cube.

    Uses ale to generate an ISD (Instrument Support Data) from the cube,
    then iterates usgscsm plugins to construct a CSM model. This replaces
    the knoten.csm.create_csm() call to avoid the knoten dependency
    (which has a broken csmapi conda dep).

    The same ISD is parsed in two ways: once as a ``csmapi.Isd`` object to
    build the CSM model, once as a Python dict to populate a
    :class:`TargetBody` with the body's ellipsoid and NAIF ID. Returning
    both together makes csm2map body-agnostic — no downstream code has to
    re-open the cube or query SPICE for radii.

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
    model : csmapi.RasterGM
        A CSM raster ground-to-image sensor model.
    body : TargetBody
        Target body ellipsoid + NAIF ID, parsed from the same ISD that
        built the model. Use this to feed ``DemRadiusSampler.fallback_radius``
        and to build body-agnostic projection strings.

    Raises
    ------
    RuntimeError
        If ale cannot generate an ISD or usgscsm cannot build the model.
    """
    import ale
    import csmapi

    from isistools.io.cubes import read_label

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
    else:
        isd_str = isd_path.read_text()

    # Parse the ISD as a dict so we can extract body info without reading
    # the file a second time. csmapi.Isd wants a filename, not a string,
    # so we can't skip the write above — but we can skip re-reading the
    # JSON file from disk by parsing the in-memory string here.
    isd_dict = json.loads(isd_str)

    # Pull the human-readable target name from the cube's Instrument
    # group when available; fall back to "BODY<naif_id>" inside TargetBody.
    target_name: str | None = None
    try:
        label = read_label(cube_path)
        inst = label["IsisCube"]["Instrument"]
        raw_target = inst.get("TargetName")
        if raw_target is not None:
            target_name = str(raw_target).strip()
    except (KeyError, FileNotFoundError, OSError):
        target_name = None

    body = TargetBody.from_isd(isd_dict, target_name=target_name)

    isd = csmapi.Isd(cube_str)

    # Iterate registered CSM plugins and try to construct a model
    for plugin in csmapi.Plugin.getList():
        for model_index in range(plugin.getNumModels()):
            model_name = plugin.getModelName(model_index)
            warnings = csmapi.WarningList()
            if plugin.canModelBeConstructedFromISD(isd, model_name, warnings):
                model = plugin.constructModelFromISD(isd, model_name)
                if isinstance(model, csmapi.RasterGM):
                    return model, body

    msg = f"No CSM plugin could construct a model from {cube_path}"
    raise RuntimeError(msg)


def load_camera_from_label(
    label_path: str | Path,
    refresh_isd: bool = True,
    use_web: bool = True,
) -> tuple["csmapi.RasterGM", "TargetBody"]:
    """Load a CSM camera model from any PDS label (EDR or cube).

    Uses ALE's NaifSpice driver to generate an ISD directly from the
    label.  By default uses SpiceQL's web service for kernel data, so
    no local SPICE kernels are needed.  No spiceinit required — this
    works with raw PDS EDR files.

    Parameters
    ----------
    label_path : path-like
        Path to a PDS3 label file (.IMG, .lbl) or ISIS cube (.cub).
        ALE will auto-detect the instrument and select the appropriate
        driver.
    refresh_isd : bool
        If True (default), regenerate the ISD even if a cached JSON exists.
    use_web : bool
        If True (default), use SpiceQL web service for SPICE data.
        Set to False to use locally installed NAIF kernels.

    Returns
    -------
    model : csmapi.RasterGM
        CSM sensor model.
    body : TargetBody
        Target body ellipsoid + NAIF ID.

    Raises
    ------
    RuntimeError
        If ALE cannot generate an ISD or no CSM plugin can build a model.
    """
    import ale
    import csmapi

    _ensure_csm_plugin_loaded()

    label_path = Path(label_path)
    label_str = str(label_path)

    isd_path = label_path.with_suffix(".json")
    if refresh_isd or not isd_path.exists():
        props = {"web": True} if use_web else {}
        isd_str = ale.loads(label_str, props=props, only_naif_spice=True)
        isd_path.write_text(isd_str)
    else:
        isd_str = isd_path.read_text()

    isd_dict = json.loads(isd_str)

    # Extract target name from the label
    target_name: str | None = None
    try:
        import pvl

        label = pvl.load(label_str)
        # Try ISIS cube format first, then PDS3
        if "IsisCube" in label:
            target_name = str(label["IsisCube"]["Instrument"].get("TargetName", ""))
        else:
            target_name = str(label.get("TARGET_NAME", ""))
        target_name = target_name.strip() or None
    except Exception:
        target_name = None

    body = TargetBody.from_isd(isd_dict, target_name=target_name)

    isd = csmapi.Isd(label_str)

    for plugin in csmapi.Plugin.getList():
        for model_index in range(plugin.getNumModels()):
            model_name = plugin.getModelName(model_index)
            warnings = csmapi.WarningList()
            if plugin.canModelBeConstructedFromISD(isd, model_name, warnings):
                model = plugin.constructModelFromISD(isd, model_name)
                if isinstance(model, csmapi.RasterGM):
                    return model, body

    msg = f"No CSM plugin could construct a model from {label_path}"
    raise RuntimeError(msg)


def get_image_size(model: "csmapi.RasterGM") -> tuple[int, int]:
    """Return (n_lines, n_samples) from a CSM model."""
    size = model.getImageSize()
    return int(size.line), int(size.samp)


def compute_ground_sample_distance(
    model: "csmapi.RasterGM",
    body: TargetBody,
) -> float:
    """Estimate the ground sample distance (GSD) in meters/pixel.

    Evaluates the CSM model at the image center: projects two adjacent
    pixels to ground, computes the great-circle distance between them
    on the target body, and returns the average of the line-direction
    and sample-direction GSDs.

    This replicates what ISIS ``cam2map`` does when no ``PixelResolution``
    is specified in the MAP file — it picks the camera's native
    resolution so the projected output neither up-samples nor
    down-samples the input.

    Parameters
    ----------
    model : csmapi.RasterGM
        CSM sensor model.
    body : TargetBody
        Target body (used for mean radius in the distance calculation).

    Returns
    -------
    float
        Ground sample distance in meters/pixel.
    """
    import csmapi

    size = model.getImageSize()
    center_line = size.line / 2.0
    center_samp = size.samp / 2.0

    # Project three points: center, center+1 in line direction,
    # center+1 in sample direction.
    ip_center = csmapi.ImageCoord(center_line, center_samp)
    ip_line = csmapi.ImageCoord(center_line + 1.0, center_samp)
    ip_samp = csmapi.ImageCoord(center_line, center_samp + 1.0)

    gp_center = model.imageToGround(ip_center, 0.0)
    gp_line = model.imageToGround(ip_line, 0.0)
    gp_samp = model.imageToGround(ip_samp, 0.0)

    def _ecef_distance(a: "csmapi.EcefCoord", b: "csmapi.EcefCoord") -> float:
        """Euclidean distance between two ECEF points, then project onto
        the sphere to get a surface arc distance. For sub-km separations
        (one pixel) the chord ≈ the arc to ~1e-8 relative error, so
        Euclidean is fine."""
        dx = a.x - b.x
        dy = a.y - b.y
        dz = a.z - b.z
        return (dx * dx + dy * dy + dz * dz) ** 0.5

    gsd_line = _ecef_distance(gp_center, gp_line)
    gsd_samp = _ecef_distance(gp_center, gp_samp)

    return (gsd_line + gsd_samp) / 2.0


def ground_to_image_batch(
    model: "csmapi.RasterGM",
    lats: np.ndarray,
    lons: np.ndarray,
    radii: np.ndarray,
    workers: int | None = None,
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
    workers : int, optional
        Number of worker threads for the per-point CSM loop. CSM calls
        release the Python GIL (they're C++), so threading gives close
        to linear speedup up to the number of CPU cores. ``None``
        (default) picks ``os.cpu_count()``. Pass ``1`` to force a
        single-threaded loop (useful for debugging and reproducibility).

    Returns
    -------
    lines, samples : ndarray
        Input image coordinates.  Invalid points get NaN.
    """
    import csmapi

    flat = lats.ravel()
    flon = lons.ravel()
    frad = radii.ravel()

    # Convert planetocentric lat/lon/radius to body-fixed ECEF (X, Y, Z).
    # This is vectorized; the slow part is the per-point CSM call below.
    cos_lat = np.cos(flat)
    sin_lat = np.sin(flat)
    cos_lon = np.cos(flon)
    sin_lon = np.sin(flon)

    x = frad * cos_lat * cos_lon
    y = frad * cos_lat * sin_lon
    z = frad * sin_lat

    n = len(flat)
    # float32 is enough for image coordinates: the line/sample values
    # sit in 0..n_lines and 0..n_samples (CTX ≤ 52 k lines × 5 k
    # samples, HiRISE ≤ 80 k × 5 k).  float32 represents integers up
    # to 16_777_216 exactly; for CTX this gives >0.001-pixel precision
    # anywhere in the frame, far below the CSM model's own pointing
    # uncertainty (sub-pixel but many multiples of 10⁻³ px).  Using
    # float32 also keeps the dtype consistent end-to-end: the
    # downstream ``_bilinear_upsample_pair`` upsamples in float32
    # already, so staying in float32 here avoids an implicit cast and
    # halves the coarse-grid memory footprint.
    out_lines = np.full(n, np.nan, dtype=np.float32)
    out_samps = np.full(n, np.nan, dtype=np.float32)

    if workers is None:
        workers = os.cpu_count() or 1

    def _process_range(i0: int, i1: int) -> None:
        """Evaluate CSM groundToImage for indices [i0, i1)."""
        for i in range(i0, i1):
            if np.isnan(flat[i]):
                continue
            try:
                gp = csmapi.EcefCoord(x[i], y[i], z[i])
                ip = model.groundToImage(gp)
                out_lines[i] = ip.line
                out_samps[i] = ip.samp
            except Exception:
                pass  # leave as NaN

    if workers <= 1 or n < 1024:
        # Single-threaded path for small batches and debugging
        _process_range(0, n)
    else:
        # Thread the per-point loop. CSM's groundToImage releases the GIL
        # in its C++ call, so threading gives near-linear speedup up to
        # the CPU count. Split the range into ``workers`` contiguous
        # chunks rather than per-item submission to minimize scheduling
        # overhead.
        from concurrent.futures import ThreadPoolExecutor

        chunk = (n + workers - 1) // workers
        ranges = [(i, min(i + chunk, n)) for i in range(0, n, chunk)]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(lambda r: _process_range(*r), ranges))

    return (
        out_lines.reshape(lats.shape),
        out_samps.reshape(lats.shape),
    )
