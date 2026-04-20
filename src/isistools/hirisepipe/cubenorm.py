"""Column normalization — Python replacement for ISIS cubenorm.

Normalizes column-to-column variations by dividing each column by
its median (or mean), preserving the global image average.
"""

from __future__ import annotations

import numpy as np


def cubenorm(
    image: np.ndarray,
    *,
    mode: str = "DIVIDE",
    normalizer: str = "MEDIAN",
    preserve: bool = True,
) -> np.ndarray:
    """Normalize column-to-column variations in a HiRISE image.

    Parameters
    ----------
    image : np.ndarray
        2D float32 image.
    mode : str
        "DIVIDE" or "SUBTRACT".
    normalizer : str
        "MEDIAN" or "AVERAGE" — per-column statistic to normalize by.
    preserve : bool
        If True, scale output to preserve the global mean.

    Returns
    -------
    np.ndarray
        Normalized image.
    """
    out = image.copy()
    n_lines, n_samps = image.shape

    # Compute per-column statistics
    if normalizer == "MEDIAN":
        col_stats = np.nanmedian(image, axis=0)
    else:
        col_stats = np.nanmean(image, axis=0)

    # Global statistic for preservation
    if preserve:
        global_stat = np.nanmedian(col_stats) if normalizer == "MEDIAN" else np.nanmean(col_stats)

    valid_cols = np.isfinite(col_stats) & (col_stats != 0)

    if mode == "DIVIDE":
        for j in range(n_samps):
            if valid_cols[j]:
                valid = np.isfinite(out[:, j])
                if preserve:
                    out[valid, j] = (out[valid, j] / col_stats[j]) * global_stat
                else:
                    out[valid, j] = out[valid, j] / col_stats[j]
    else:
        for j in range(n_samps):
            if valid_cols[j]:
                valid = np.isfinite(out[:, j])
                if preserve:
                    out[valid, j] = out[valid, j] - col_stats[j] + global_stat
                else:
                    out[valid, j] = out[valid, j] - col_stats[j]

    return out
