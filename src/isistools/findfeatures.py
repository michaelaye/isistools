"""Feature matching for image coregistration — Python replacement for ISIS findfeatures.

Uses OpenCV's AKAZE (or ORB/SIFT) detector and BFMatcher with Lowe's ratio
test to find tie points between image pairs.  Outputs a pandas DataFrame
compatible with plio's ``to_isis()`` for writing ISIS binary control networks.

Usage::

    from isistools.findfeatures import find_features, match_pair

    # From numpy arrays
    matches = match_pair(image1, image2)

    # From files, with control network output
    cnet = find_features(
        from_path="image1.cub",
        match_path="image2.cub",
        onet_path="output.net",
    )
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np
import pandas as pd


@dataclass
class MatchResult:
    """Result of feature matching between two images."""

    # Matched keypoint coordinates (sample, line) in each image
    from_samples: np.ndarray
    from_lines: np.ndarray
    match_samples: np.ndarray
    match_lines: np.ndarray

    # Match quality metrics
    distances: np.ndarray  # descriptor distances for each match
    n_keypoints_from: int  # total keypoints detected in FROM image
    n_keypoints_match: int  # total keypoints detected in MATCH image
    n_matches_raw: int  # matches before ratio test
    n_matches_good: int  # matches after ratio test

    @property
    def n_points(self) -> int:
        return len(self.from_samples)


def _normalize_image(image: np.ndarray) -> np.ndarray:
    """Normalize a float image to uint8 for OpenCV feature detection."""
    if image.dtype == np.uint8:
        return image

    # Handle NaN
    valid = np.isfinite(image)
    if not valid.any():
        return np.zeros(image.shape, dtype=np.uint8)

    # Percentile stretch to uint8
    vmin = np.nanpercentile(image, 0.5)
    vmax = np.nanpercentile(image, 99.5)
    if vmax == vmin:
        return np.zeros(image.shape, dtype=np.uint8)

    stretched = np.clip((image - vmin) / (vmax - vmin) * 255, 0, 255)
    stretched = np.nan_to_num(stretched, nan=0).astype(np.uint8)
    return stretched


def match_pair(
    image_from: np.ndarray,
    image_match: np.ndarray,
    *,
    algorithm: str = "AKAZE",
    ratio: float = 0.65,
    max_points: int = 0,
) -> MatchResult:
    """Find matching features between two images.

    Parameters
    ----------
    image_from : np.ndarray
        First image (2D array, any numeric dtype).
    image_match : np.ndarray
        Second image (2D array, any numeric dtype).
    algorithm : str
        Feature detector: "AKAZE" (default), "ORB", or "SIFT".
    ratio : float
        Lowe's ratio test threshold (default 0.65, matching ISIS default).
        Lower = stricter matching, fewer but more reliable matches.
    max_points : int
        Maximum number of matches to return (0 = unlimited).
        If set, the best matches (by distance) are kept.

    Returns
    -------
    MatchResult
        Matched keypoint coordinates and quality metrics.
    """
    # Normalize to uint8 for OpenCV
    img1 = _normalize_image(image_from)
    img2 = _normalize_image(image_match)

    # Create feature detector
    if algorithm.upper() == "AKAZE":
        detector = cv2.AKAZE_create()
        norm_type = cv2.NORM_HAMMING
    elif algorithm.upper() == "ORB":
        detector = cv2.ORB_create(nfeatures=50000)
        norm_type = cv2.NORM_HAMMING
    elif algorithm.upper() == "SIFT":
        detector = cv2.SIFT_create()
        norm_type = cv2.NORM_L2
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}. Use AKAZE, ORB, or SIFT.")

    # Detect and compute descriptors
    kp1, desc1 = detector.detectAndCompute(img1, None)
    kp2, desc2 = detector.detectAndCompute(img2, None)

    if desc1 is None or desc2 is None or len(kp1) == 0 or len(kp2) == 0:
        return MatchResult(
            from_samples=np.array([]),
            from_lines=np.array([]),
            match_samples=np.array([]),
            match_lines=np.array([]),
            distances=np.array([]),
            n_keypoints_from=len(kp1) if kp1 else 0,
            n_keypoints_match=len(kp2) if kp2 else 0,
            n_matches_raw=0,
            n_matches_good=0,
        )

    # Match descriptors with kNN (k=2 for ratio test)
    matcher = cv2.BFMatcher(norm_type)
    raw_matches = matcher.knnMatch(desc1, desc2, k=2)
    n_raw = len(raw_matches)

    # Apply Lowe's ratio test
    good_matches = []
    for match_pair_result in raw_matches:
        if len(match_pair_result) < 2:
            continue
        m, n = match_pair_result
        if m.distance < ratio * n.distance:
            good_matches.append(m)

    # Sort by distance (best first)
    good_matches.sort(key=lambda m: m.distance)

    # Limit if requested
    if max_points > 0 and len(good_matches) > max_points:
        good_matches = good_matches[:max_points]

    # Extract coordinates
    from_pts = np.array([kp1[m.queryIdx].pt for m in good_matches])
    match_pts = np.array([kp2[m.trainIdx].pt for m in good_matches])
    distances = np.array([m.distance for m in good_matches])

    if len(from_pts) == 0:
        from_pts = np.empty((0, 2))
        match_pts = np.empty((0, 2))

    return MatchResult(
        from_samples=from_pts[:, 0] if len(from_pts) > 0 else np.array([]),
        from_lines=from_pts[:, 1] if len(from_pts) > 0 else np.array([]),
        match_samples=match_pts[:, 0] if len(match_pts) > 0 else np.array([]),
        match_lines=match_pts[:, 1] if len(match_pts) > 0 else np.array([]),
        distances=distances,
        n_keypoints_from=len(kp1),
        n_keypoints_match=len(kp2),
        n_matches_raw=n_raw,
        n_matches_good=len(good_matches),
    )


def _image_to_ground(model, sample: float, line: float) -> tuple[float, float, float]:
    """Convert image coordinates to body-fixed XYZ using a CSM model.

    Parameters
    ----------
    model : csmapi.RasterGM
        CSM sensor model.
    sample, line : float
        Image coordinates (0-based, pixel center at 0.0).

    Returns
    -------
    tuple of (x, y, z) in meters (body-fixed ECEF).
    """
    import csmapi

    ic = csmapi.ImageCoord(line, sample)
    try:
        gc = model.imageToGround(ic, 0.0)  # height = 0 (on ellipsoid)
        return (gc.x, gc.y, gc.z)
    except Exception:
        return (0.0, 0.0, 0.0)


def _compute_ground_points(
    model,
    samples: np.ndarray,
    lines: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute body-fixed XYZ for arrays of image coordinates.

    Returns (x_array, y_array, z_array) in meters.
    """
    n = len(samples)
    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    for i in range(n):
        x[i], y[i], z[i] = _image_to_ground(model, samples[i], lines[i])
    return x, y, z


def matches_to_cnet(
    result: MatchResult,
    serial_from: str,
    serial_match: str,
    *,
    model_from=None,
    point_id_prefix: str = "feature_",
    point_id_start: int = 1,
) -> pd.DataFrame:
    """Convert match results to a plio-compatible control network DataFrame.

    Parameters
    ----------
    result : MatchResult
        Output from ``match_pair()``.
    serial_from : str
        ISIS serial number of the FROM image.
    serial_match : str
        ISIS serial number of the MATCH image.
    model_from : csmapi.RasterGM, optional
        CSM camera model for the FROM image.  If provided, computes
        body-fixed XYZ ground coordinates (aprioriX/Y/Z) for each
        point using the reference measure's image coordinates.
        Without this, ground coordinates are set to (0, 0, 0) and
        jigsaw will not be able to use the network.
    point_id_prefix : str
        Prefix for point IDs (e.g. "feature_" → "feature_00001").
    point_id_start : int
        Starting index for point IDs.

    Returns
    -------
    pd.DataFrame
        Control network DataFrame with columns compatible with
        ``plio.io.io_controlnetwork.to_isis()``.
    """
    # Compute ground coordinates from the FROM (reference) image if model provided
    if model_from is not None:
        gx, gy, gz = _compute_ground_points(model_from, result.from_samples, result.from_lines)
    else:
        gx = np.zeros(result.n_points)
        gy = np.zeros(result.n_points)
        gz = np.zeros(result.n_points)

    rows = []
    for i in range(result.n_points):
        pid = f"{point_id_prefix}{point_id_start + i:05d}"

        # Shared point-level fields (ground coordinates from reference measure)
        point_fields = {
            "id": pid,
            "pointType": 2,  # Free
            "aprioriX": gx[i],
            "aprioriY": gy[i],
            "aprioriZ": gz[i],
            "adjustedX": gx[i],
            "adjustedY": gy[i],
            "adjustedZ": gz[i],
            "aprioriSurfPointSource": 1 if model_from is not None else 0,
            "referenceIndex": 0,
        }

        # FROM measure (reference)
        rows.append(
            {
                **point_fields,
                "serialnumber": serial_from,
                "measureType": 3,  # Reference
                "sample": float(result.from_samples[i]),
                "line": float(result.from_lines[i]),
                "sampleResidual": 0.0,
                "lineResidual": 0.0,
            }
        )

        # MATCH measure (Candidate)
        rows.append(
            {
                **point_fields,
                "serialnumber": serial_match,
                "measureType": 0,  # Candidate
                "sample": float(result.match_samples[i]),
                "line": float(result.match_lines[i]),
                "sampleResidual": 0.0,
                "lineResidual": 0.0,
            }
        )

    return pd.DataFrame(rows)


def _read_image(path: str | Path) -> np.ndarray:
    """Read an image from a file (ISIS cube, GeoTIFF, or standard image)."""
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in (".cub",):
        import warnings

        import rioxarray  # noqa: F401
        import xarray as xr

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="Dataset has no geotransform")
            da = xr.open_dataarray(str(path), engine="rasterio")
        return da.values.squeeze()
    elif suffix in (".tif", ".tiff"):
        import rasterio

        with rasterio.open(str(path)) as src:
            return src.read(1)
    else:
        # Try OpenCV for standard image formats
        img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        if img is None:
            raise ValueError(f"Cannot read image: {path}")
        return img


def _get_serial_number(cube_path: str | Path) -> str:
    """Extract or construct an ISIS serial number from a cube label.

    For MRO instruments, the serial number format is:
    MRO/{INSTRUMENT}/{SpacecraftClockCount}

    For other instruments, falls back to the product ID.
    """

    from isistools.io.cubes import read_label

    try:
        label = read_label(cube_path)
        inst = label["IsisCube"]["Instrument"]

        spacecraft = str(inst.get("SpacecraftName", "")).replace(" ", "_").upper()
        instrument_id = str(inst.get("InstrumentId", ""))
        sclk = str(inst.get("SpacecraftClockCount", inst.get("SpacecraftClockStartCount", "")))

        if "MARS_RECONNAISSANCE" in spacecraft:
            return f"MRO/{instrument_id}/{sclk}"
        else:
            return f"{spacecraft}/{instrument_id}/{sclk}"
    except Exception:
        # Fallback: use filename
        return str(Path(cube_path).stem)


def find_features(
    from_path: str | Path,
    match_path: str | Path,
    *,
    onet_path: str | Path | None = None,
    model_from=None,
    from_cube: str | Path | None = None,
    algorithm: str = "AKAZE",
    ratio: float = 0.65,
    max_points: int = 0,
    target_name: str = "Mars",
    network_id: str = "Features",
) -> pd.DataFrame:
    """Find feature matches between two images and optionally write a control network.

    This is the high-level function replacing ISIS ``findfeatures``.

    Parameters
    ----------
    from_path : path-like
        Path to the FROM image (ISIS cube, GeoTIFF, or standard image).
    match_path : path-like
        Path to the MATCH image.
    onet_path : path-like, optional
        If provided, write the control network as an ISIS binary .net file.
    model_from : csmapi.RasterGM, optional
        CSM camera model for the FROM image.  If provided, ground
        coordinates (aprioriX/Y/Z) are computed for each tie point.
    from_cube : path-like, optional
        Path to the FROM image's spiceinit'd ISIS cube.  If ``model_from``
        is not provided but ``from_cube`` is, the camera model is loaded
        from the cube via ``csm2map.camera.load_camera()``.
    algorithm : str
        Feature detector: "AKAZE" (default), "ORB", or "SIFT".
    ratio : float
        Lowe's ratio test threshold (default 0.65).
    max_points : int
        Maximum number of matches (0 = unlimited).
    target_name : str
        Target body name for the control network (default "Mars").
    network_id : str
        Network ID string for the control network.

    Returns
    -------
    pd.DataFrame
        Control network DataFrame with matched tie points.
    """
    from_path = Path(from_path)
    match_path = Path(match_path)

    # Load camera model if cube path provided
    if model_from is None and from_cube is not None:
        from isistools.csm2map.camera import load_camera

        model_from, _body = load_camera(from_cube)

    # Read images
    img_from = _read_image(from_path)
    img_match = _read_image(match_path)

    # Get serial numbers
    serial_from = _get_serial_number(from_path)
    serial_match = _get_serial_number(match_path)

    # Match features
    result = match_pair(
        img_from,
        img_match,
        algorithm=algorithm,
        ratio=ratio,
        max_points=max_points,
    )

    # Convert to control network DataFrame (with ground coords if model available)
    cnet = matches_to_cnet(result, serial_from, serial_match, model_from=model_from)

    # Write if requested
    if onet_path is not None and len(cnet) > 0:
        from plio.io.io_controlnetwork import to_isis

        to_isis(
            cnet,
            str(onet_path),
            targetname=target_name,
            networkid=network_id,
        )

    return cnet
