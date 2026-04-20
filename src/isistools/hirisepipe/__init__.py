"""hirisepipe — Python HiRISE RED mosaic pipeline.

A pure-Python replacement for the ISIS3 HiRISE RED processing chain::

    hi2isis → hical → histitch → cubenorm → cam2map → equalizer → automos

Quick start
-----------
Python::

    from isistools.hirisepipe import hical, stitch_channels, cubenorm
    cal = hical(raw_cube_path, units="DN")

Public API
----------
hical : function
    Radiometric calibration (replaces ISIS hical).
stitch_channels : function
    Channel stitching with balance (replaces ISIS histitch).
cubenorm : function
    Column normalization (replaces ISIS cubenorm).
process_ccd : function
    Full single-CCD pipeline: hical → histitch → cubenorm.
"""

from isistools.hirisepipe.cubenorm import cubenorm
from isistools.hirisepipe.hical import hical, hical_from_edr
from isistools.hirisepipe.ingest import HiRISEEDR, ingest_hirise_edr
from isistools.hirisepipe.pipeline import calibrate_all, calibrate_ccd, create_red_mosaic
from isistools.hirisepipe.stitch import stitch_channels

__all__ = [
    "HiRISEEDR",
    "create_red_mosaic",
    "cubenorm",
    "hical",
    "hical_from_edr",
    "ingest_hirise_edr",
    "calibrate_all",
    "calibrate_ccd",
    "stitch_channels",
]
