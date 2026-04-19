"""ctxpipe — Python CTX calibration and projection pipeline.

A pure-Python replacement for the ISIS3 CTX processing chain::

    mroctx2isis → ctxcal → ctxevenodd → cam2map

Products
--------
- **Level 0** (raw): ``ingest_ctx_edr()`` — SQROOT-decoded numpy array
- **Level 1** (calibrated): ``ctx_calibrate()`` — dark/flat/evenodd corrected
- **Level 2** (projected): ``ctx_project()`` or ``ctx_edr_to_map()`` — GeoTIFF

Quick start
-----------
Calibrate only::

    from isistools.ctxpipe import ctx_calibrate
    cal, meta = ctx_calibrate("input.IMG")

Calibrate without even/odd correction::

    cal, meta = ctx_calibrate("input.IMG", evenodd=False)

Full pipeline to map-projected GeoTIFF::

    from isistools.ctxpipe import ctx_edr_to_map
    ctx_edr_to_map("input.IMG", "input.cub", "projected.tif")
"""

from isistools.ctxpipe.calibrate import calibrate
from isistools.ctxpipe.evenodd import correct_evenodd
from isistools.ctxpipe.ingest import ingest_ctx_edr
from isistools.ctxpipe.pipeline import ctx_calibrate, ctx_edr_to_map, ctx_project

__all__ = [
    "calibrate",
    "correct_evenodd",
    "ctx_calibrate",
    "ctx_edr_to_map",
    "ctx_project",
    "ingest_ctx_edr",
]
