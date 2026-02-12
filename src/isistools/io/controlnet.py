"""Control network I/O using plio.

Wraps plio's control network reader/writer and adds convenience
columns for visualization (residual magnitude, point status classification).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from plio.io.io_controlnetwork import from_isis, to_isis


def _classify_point_status(row: pd.Series) -> str:
    """Classify a control net measure into a display status.

    Categories
    ----------
    - ``ignored``: point or measure is flagged as ignored
    - ``registered``: measure has been successfully registered
      (has non-zero residuals or is of type RegisteredPixel/RegisteredSubPixel)
    - ``unregistered``: candidate measure not yet registered
    """
    if row.get("pointIgnore", False) or row.get("measureIgnore", False):
        return "ignored"

    measure_type = row.get("measureType", 0)
    # plio measure types: 0=Candidate, 1=Manual, 2=RegisteredPixel,
    # 3=RegisteredSubPixel
    if measure_type >= 2:
        return "registered"

    # Check if residuals exist (non-zero means it went through jigsaw)
    res_s = row.get("residualSample", 0.0)
    res_l = row.get("residualLine", 0.0)
    if abs(res_s) > 1e-10 or abs(res_l) > 1e-10:
        return "registered"

    return "unregistered"


def load_cnet(path: str | Path) -> pd.DataFrame:
    """Load an ISIS control network file.

    Uses plio to read the binary protobuf format and adds
    convenience columns for visualization.

    Parameters
    ----------
    path : path-like
        Path to the .net control network file.

    Returns
    -------
    pd.DataFrame
        Control network with original plio columns plus:
        - ``residual_magnitude``: Euclidean residual
        - ``status``: 'registered', 'unregistered', or 'ignored'

    Notes
    -----
    The returned DataFrame has one row per *measure* (not per point).
    Multiple measures share the same ``pointId``. The ``serialnumber``
    column identifies which cube image the measure belongs to.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Control network not found: {path}")

    from isistools.io.cache import get_cache

    cache = get_cache()
    cache_key = f"cnet:{path}:{path.stat().st_mtime_ns}"
    cached = cache.get(cache_key)
    if cached is not None:
        return cached

    df = from_isis(str(path))

    # Normalize plio column names to isistools conventions
    df = df.rename(columns={
        "id": "pointId",
        "sampleResidual": "residualSample",
        "lineResidual": "residualLine",
    })

    # Add residual magnitude
    res_s = df.get("residualSample", pd.Series(0.0, index=df.index))
    res_l = df.get("residualLine", pd.Series(0.0, index=df.index))
    df["residual_magnitude"] = np.hypot(
        res_s.fillna(0.0).astype(float),
        res_l.fillna(0.0).astype(float),
    )

    # Classify status
    df["status"] = df.apply(_classify_point_status, axis=1)

    # Keep only columns isistools uses — plio's IsisControlNetwork carries
    # protobuf repeated-message fields that cannot be pickled.
    _keep = [
        "pointId", "serialnumber", "sample", "line",
        "residualSample", "residualLine", "residual_magnitude",
        "measureType", "pointType", "pointIgnore", "measureIgnore", "status",
        "adjustedX", "adjustedY", "adjustedLon", "adjustedLat",
        "aprioriX", "aprioriY", "aprioriLon", "aprioriLat",
    ]
    df = pd.DataFrame(df[[c for c in _keep if c in df.columns]])

    cache.set(cache_key, df)
    return df


def save_cnet(df: pd.DataFrame, path: str | Path) -> None:
    """Write a control network DataFrame back to ISIS binary format.

    Parameters
    ----------
    df : pd.DataFrame
        Control network DataFrame (as returned by :func:`load_cnet`).
    path : path-like
        Output path for the .net file.
    """
    to_isis(df, str(path))


def cnet_summary(df: pd.DataFrame) -> dict:
    """Compute summary statistics for a control network.

    Returns
    -------
    dict
        Keys: n_points, n_measures, n_images, n_registered,
        n_unregistered, n_ignored, mean_residual, max_residual.
    """
    return {
        "n_points": df["pointId"].nunique() if "pointId" in df.columns else 0,
        "n_measures": len(df),
        "n_images": df["serialnumber"].nunique() if "serialnumber" in df.columns else 0,
        "n_registered": (df["status"] == "registered").sum(),
        "n_unregistered": (df["status"] == "unregistered").sum(),
        "n_ignored": (df["status"] == "ignored").sum(),
        "mean_residual": df["residual_magnitude"].mean(),
        "max_residual": df["residual_magnitude"].max(),
    }
