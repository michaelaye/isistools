"""CTX radiometric calibration — Python replacement for ISIS ctxcal.

Applies dark current subtraction, flat-field correction, and optional
I/F conversion.  Follows the ISIS ctxcal algorithm exactly.
"""

from __future__ import annotations

import os
import warnings
from pathlib import Path

import numpy as np

from isistools.ctxpipe.ingest import CTXMetadata
from isistools.special_pixels import VALID_MIN


def _is_special(pixel: np.ndarray) -> np.ndarray:
    """Test if pixels are ISIS special sentinel values (not NaN)."""
    # Compare in float64 to avoid float32 overflow warnings
    return pixel.astype(np.float64) <= VALID_MIN


# Default calibration data directory
_DEFAULT_ISISDATA = Path(
    os.environ.get("ISISDATA", str(Path.home() / "Dropbox" / "data" / "isisdata"))
)

# CTX calibration constants (from ISIS ctxcal source)
# w0: instrument sensitivity in DN/msec at 1 AU equivalent
_W0 = 3660.5
# Mars perihelion distance in km (used for I/F normalization)
_PERIHELION_KM = 2.07e8


def _find_calibration_file(
    pattern: str,
    calibration_dir: Path | None = None,
) -> Path:
    """Find the highest-versioned calibration file matching a pattern."""
    if calibration_dir is None:
        calibration_dir = _DEFAULT_ISISDATA / "mro" / "calibration"

    candidates = sorted(calibration_dir.glob(pattern))
    if not candidates:
        raise FileNotFoundError(f"No calibration file matching '{pattern}' in {calibration_dir}")
    return candidates[-1]  # highest version


def _load_flat_field(
    calibration_dir: Path | None = None,
    flat_path: Path | None = None,
) -> np.ndarray:
    """Load the CTX flat-field calibration cube.

    Returns
    -------
    flat : np.ndarray
        1D float32 array of 5000 samples (full detector width).
    """
    if flat_path is None:
        flat_path = _find_calibration_file("ctxFlat_????.cub", calibration_dir)

    # The flat field is a 1-band, 5000-sample, 1-line ISIS cube.
    # Read it via rioxarray/GDAL.
    import rioxarray  # noqa: F401
    import xarray as xr

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Dataset has no geotransform")
        da = xr.open_dataarray(flat_path, engine="rasterio")

    flat = da.values.ravel().astype(np.float32)
    if flat.shape[0] != 5000:
        raise ValueError(f"Expected 5000-sample flat field, got {flat.shape[0]}")
    return flat


def _compute_dark_current(
    dark_pixels: np.ndarray,
    spatial_summing: int,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Compute per-line dark current from prefix dark pixels.

    Parameters
    ----------
    dark_pixels : np.ndarray
        Shape (n_lines, n_dark_pixels), int16.  -32768 = NULL.
    spatial_summing : int
        1 or 2.

    Returns
    -------
    dc_a : np.ndarray
        Per-line dark current for channel A (or combined if summing > 1).
        Shape (n_lines,).
    dc_b : np.ndarray or None
        Per-line dark current for channel B (None if summing > 1).
    """
    n_lines = dark_pixels.shape[0]
    valid_mask = dark_pixels != -32768  # not ISIS NULL

    if spatial_summing == 1:
        # Channels alternate: A, B, A, B, ... (odd/even columns)
        # Column 0 = A, column 1 = B, column 2 = A, ...
        a_mask = np.zeros(dark_pixels.shape[1], dtype=bool)
        b_mask = np.zeros(dark_pixels.shape[1], dtype=bool)
        a_mask[0::2] = True
        b_mask[1::2] = True

        dc_a = np.full(n_lines, np.nan, dtype=np.float64)
        dc_b = np.full(n_lines, np.nan, dtype=np.float64)

        for i in range(n_lines):
            a_valid = valid_mask[i] & a_mask
            b_valid = valid_mask[i] & b_mask

            if a_valid.any():
                dc_a[i] = dark_pixels[i, a_valid].astype(np.float64).mean()
            if b_valid.any():
                dc_b[i] = dark_pixels[i, b_valid].astype(np.float64).mean()

        return dc_a, dc_b
    else:
        # Summing > 1: channels are mixed, use single average
        dc = np.full(n_lines, np.nan, dtype=np.float64)
        for i in range(n_lines):
            line_valid = valid_mask[i]
            if line_valid.any():
                dc[i] = dark_pixels[i, line_valid].astype(np.float64).mean()
        return dc, None


def calibrate(
    image: np.ndarray,
    metadata: CTXMetadata,
    *,
    iof: bool = False,
    sun_distance_km: float | None = None,
    calibration_dir: Path | None = None,
    flat_path: Path | None = None,
) -> np.ndarray:
    """Apply CTX radiometric calibration.

    Parameters
    ----------
    image : np.ndarray
        2D float32 array from ``ingest_ctx_edr``.
    metadata : CTXMetadata
        Metadata from ``ingest_ctx_edr``, including dark_pixels.
    iof : bool
        If True, convert to I/F.  Requires ``sun_distance_km``.
        If False (default), output is DN/ms.
    sun_distance_km : float, optional
        Sun-to-target distance in km.  Required if ``iof=True``.
        Can be computed from SPICE (e.g., via spiceypy).
    calibration_dir : Path, optional
        Directory containing CTX calibration files.
        Defaults to ``~/Dropbox/data/isisdata/mro/calibration``.
    flat_path : Path, optional
        Explicit path to the flat-field cube.

    Returns
    -------
    calibrated : np.ndarray
        2D float32 array of calibrated values.
    """
    if iof and sun_distance_km is None:
        raise ValueError("sun_distance_km is required when iof=True")

    exposure = metadata.line_exposure_duration  # ms

    # Load flat field
    flat_full = _load_flat_field(calibration_dir, flat_path)

    # Compute dark current per line
    dc_a, dc_b = _compute_dark_current(metadata.dark_pixels, metadata.spatial_summing)

    # Determine flat-field offset.
    # ISIS ctxcal: if firstSamp > 0, firstSamp -= 38
    first_samp = metadata.sample_first_pixel
    if first_samp > 0:
        first_samp -= 38

    n_lines, n_samps = image.shape

    if metadata.spatial_summing == 1:
        # Build flat field array for this image width
        flat_slice = flat_full[first_samp : first_samp + n_samps]
        if flat_slice.shape[0] < n_samps:
            raise ValueError(
                f"Flat field too short: need {n_samps} samples starting "
                f"at {first_samp}, flat has {flat_full.shape[0]}"
            )

        # Build 2D dark current array: alternating A/B channels per line
        # dark_2d[i, even_col] = dc_a[i], dark_2d[i, odd_col] = dc_b[i]
        dark_2d = np.empty((n_lines, n_samps), dtype=np.float32)
        dark_2d[:, 0::2] = dc_a[:, np.newaxis]
        if dc_b is not None:
            dark_2d[:, 1::2] = dc_b[:, np.newaxis]
        else:
            dark_2d[:, 1::2] = dc_a[:, np.newaxis]

        # Vectorized calibration: (image - dark) / (exposure * flat)
        valid = np.isfinite(image)
        out = np.where(
            valid,
            (image - dark_2d) / (exposure * flat_slice[np.newaxis, :]),
            image,
        ).astype(np.float32)
    else:
        # Summing mode 2: average adjacent flat pixels, single dark
        flat_summed = (
            flat_full[first_samp : first_samp + n_samps * 2 : 2]
            + flat_full[first_samp + 1 : first_samp + n_samps * 2 : 2]
        ) / 2.0

        # Vectorized: dc_a is per-line, broadcast over samples
        valid = np.isfinite(image)
        out = np.where(
            valid,
            (image - dc_a[:, np.newaxis]) / (exposure * flat_summed[np.newaxis, :]),
            image,
        ).astype(np.float32)

    # I/F conversion
    if iof:
        w1 = _W0 * (_PERIHELION_KM**2) / (sun_distance_km**2)
        # ISIS formula: ((DN - dark) / flat) * (1 / (exposure * w1))
        # We already have (DN - dark) / (exposure * flat), so divide by w1
        out = out / w1

    return out
