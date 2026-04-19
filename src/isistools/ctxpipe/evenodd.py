"""CTX even/odd column correction — Python replacement for ISIS ctxevenodd.

Corrects the systematic striping pattern caused by alternating column
readout in the CTX detector.  Only applies to images with spatial
summing = 1.
"""

from __future__ import annotations

import numpy as np


def correct_evenodd(
    image: np.ndarray,
    spatial_summing: int = 1,
) -> np.ndarray:
    """Apply even/odd column correction to a CTX image.

    Parameters
    ----------
    image : np.ndarray
        2D float32 array (calibrated CTX image).
    spatial_summing : int
        Spatial summing mode.  Correction is only applied for
        summing = 1.  For summing > 1, the image is returned unchanged.

    Returns
    -------
    corrected : np.ndarray
        2D float32 array with even/odd correction applied.

    Raises
    ------
    ValueError
        If the image has no valid pixels in odd or even columns.
    """
    if spatial_summing != 1:
        return image.copy()

    # Compute mean of odd and even columns using views (no copy)
    odd_cols = image[:, 0::2]
    even_cols = image[:, 1::2]

    odd_mean = float(np.nanmean(odd_cols))
    even_mean = float(np.nanmean(even_cols))

    if not np.isfinite(odd_mean):
        raise ValueError("No valid pixels in odd columns")
    if not np.isfinite(even_mean):
        raise ValueError("No valid pixels in even columns")

    # Correction offset = half the difference between even and odd means
    correction = np.float32((even_mean - odd_mean) / 2.0)

    # Apply in-place on column slices (no full-image mask allocation)
    out = image.copy()
    odd_view = out[:, 0::2]
    even_view = out[:, 1::2]

    # Only correct finite pixels
    odd_valid = np.isfinite(odd_view)
    odd_view[odd_valid] += correction

    even_valid = np.isfinite(even_view)
    even_view[even_valid] -= correction

    return out
