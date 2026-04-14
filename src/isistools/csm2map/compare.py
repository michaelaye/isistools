"""Compare csm2map output against an ISIS cam2map reference.

Reports pixel-level difference statistics for validation. Usable as a
library function (e.g. from a notebook) or via the ``csm2map compare``
CLI command.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import rasterio
from rich.console import Console

from isistools.io.cubes import read_isis_cube_raw


def compare(
    isis_projected: str | Path,
    csm_projected: str | Path,
    console: Console | None = None,
) -> dict:
    """Compare an ISIS cam2map cube with a csm2map GeoTIFF.

    Parameters
    ----------
    isis_projected : path-like
        Path to the ISIS cam2map output cube.
    csm_projected : path-like
        Path to the csm2map output GeoTIFF.
    console : rich.Console, optional
        Console for output. If None, a new one is created.

    Returns
    -------
    dict
        Comparison statistics: n_both, n_isis_only, n_csm_only,
        mean_diff, median_diff, std_diff, min_diff, max_diff, and
        a dict of threshold percentages.

    Raises
    ------
    ValueError
        If the two images have different shapes (mismatched grids).
    """
    if console is None:
        console = Console()

    console.print("[bold]Loading ISIS projected cube[/bold]")
    isis_data, _ = read_isis_cube_raw(isis_projected)

    console.print("[bold]Loading CSM projected GeoTIFF[/bold]")
    with rasterio.open(str(csm_projected)) as src:
        csm_data = src.read(1).astype(np.float32)

    if isis_data.shape != csm_data.shape:
        msg = (
            f"Shape mismatch: ISIS {isis_data.shape} vs CSM {csm_data.shape}. "
            f"Comparison requires matching grid parameters."
        )
        raise ValueError(msg)

    # Compare only where both have valid data
    isis_valid = np.isfinite(isis_data) & (isis_data != 0)
    csm_valid = np.isfinite(csm_data) & (csm_data != 0)
    both_valid = isis_valid & csm_valid

    n_both = int(np.sum(both_valid))
    n_isis_only = int(np.sum(isis_valid & ~csm_valid))
    n_csm_only = int(np.sum(csm_valid & ~isis_valid))

    console.print(f"\n  Both valid: {n_both:,}")
    console.print(f"  ISIS-only:  {n_isis_only:,}")
    console.print(f"  CSM-only:   {n_csm_only:,}")

    result: dict = {
        "n_both": n_both,
        "n_isis_only": n_isis_only,
        "n_csm_only": n_csm_only,
    }

    if n_both == 0:
        console.print("[red]No overlapping valid pixels![/red]")
        return result

    diff = csm_data[both_valid] - isis_data[both_valid]
    result.update(
        {
            "mean_diff": float(np.mean(diff)),
            "median_diff": float(np.median(diff)),
            "std_diff": float(np.std(diff)),
            "min_diff": float(np.min(diff)),
            "max_diff": float(np.max(diff)),
        }
    )

    console.print("\n  [bold]Difference statistics (CSM - ISIS):[/bold]")
    console.print(f"  Mean:   {result['mean_diff']:.4f}")
    console.print(f"  Median: {result['median_diff']:.4f}")
    console.print(f"  Std:    {result['std_diff']:.4f}")
    console.print(f"  Min:    {result['min_diff']:.4f}")
    console.print(f"  Max:    {result['max_diff']:.4f}")

    thresholds = {}
    for threshold in [0.01, 0.1, 1.0, 5.0]:
        pct = 100 * np.sum(np.abs(diff) < threshold) / n_both
        thresholds[threshold] = float(pct)
        console.print(f"  |diff| < {threshold}: {pct:.1f}%")

    result["thresholds"] = thresholds
    return result
