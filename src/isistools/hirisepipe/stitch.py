"""HiRISE channel stitching — Python replacement for ISIS histitch.

Stitches two CCD channels (0 and 1) into a single image,
optionally balancing their radiometric levels at the seam.
"""

from __future__ import annotations

import numpy as np


def stitch_channels(
    channel0: np.ndarray,
    channel1: np.ndarray,
    *,
    balance: bool = True,
    truth_channel: int = 0,
    seam_size: int = 128,
    skip: int = 5,
) -> np.ndarray:
    """Stitch two HiRISE CCD channels with optional balance correction.

    Channel 1 goes on the left, channel 0 goes on the right (HiRISE
    detector layout).

    Parameters
    ----------
    channel0, channel1 : np.ndarray
        Calibrated images for channels 0 and 1.  Shape (n_lines, n_samples).
    balance : bool
        If True, apply a multiplicative gain correction so the seam edge
        averages match.
    truth_channel : int
        Channel to hold fixed (0 or 1).  The other is scaled.
    seam_size : int
        Number of pixels at the seam edge to use for balance statistics.
    skip : int
        Number of pixels to skip from the seam edge before computing stats.

    Returns
    -------
    np.ndarray
        Stitched image, shape (n_lines, n_samples_ch0 + n_samples_ch1).
    """
    n_lines = max(channel0.shape[0], channel1.shape[0])
    n_samps_0 = channel0.shape[1]
    n_samps_1 = channel1.shape[1]

    if balance:
        # Compute overlap statistics at the seam
        # Channel 0 seam: leftmost pixels (after skip)
        ch0_seam = channel0[:, skip : skip + seam_size + 1]
        # Channel 1 seam: rightmost pixels (after skip)
        ch1_seam = channel1[:, -(skip + seam_size + 1) : -skip if skip > 0 else None]

        ch0_valid = ch0_seam[np.isfinite(ch0_seam)]
        ch1_valid = ch1_seam[np.isfinite(ch1_seam)]

        if len(ch0_valid) > 0 and len(ch1_valid) > 0:
            avg0 = ch0_valid.mean()
            avg1 = ch1_valid.mean()

            if truth_channel == 0:
                coeff = avg0 / avg1 if avg1 != 0 else 1.0
                channel1_adj = channel1.copy()
                valid = np.isfinite(channel1_adj)
                channel1_adj[valid] *= coeff
            else:
                coeff = avg1 / avg0 if avg0 != 0 else 1.0
                channel0_adj = channel0.copy()
                valid = np.isfinite(channel0_adj)
                channel0_adj[valid] *= coeff
                channel0 = channel0_adj
                channel1_adj = channel1
        else:
            channel1_adj = channel1
    else:
        channel1_adj = channel1

    # Concatenate: channel 1 on left, channel 0 on right
    out = np.full((n_lines, n_samps_0 + n_samps_1), np.nan, dtype=np.float32)
    out[: channel1_adj.shape[0], :n_samps_1] = channel1_adj
    out[: channel0.shape[0], n_samps_1:] = channel0

    return out
