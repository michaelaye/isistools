"""csm2map — CSM-based map projection for planetary images.

A drop-in replacement for ISIS ``cam2map`` using the Community Sensor
Model (CSM) standard.  Produces GeoTIFF output with ISIS-compatible PVL
sidecar metadata.  Body-agnostic (works with any target ALE supports),
handles all ISIS MAP file conventions (planetographic/planetocentric,
positive-east/west, 180/360 domain), and auto-derives resolution +
bounds from the camera model when no MAP file is given.

Quick start
-----------
CLI (zero flags)::

    csm2map input.cub output.tif

Python (one call)::

    from isistools.csm2map import csm2map
    csm2map("input.cub", "output.tif")

Step-by-step Python::

    from isistools.csm2map.camera import load_camera
    from isistools.csm2map.grid import grid_from_params
    from isistools.csm2map.transform import compute_transform_coarse
    from isistools.csm2map.resample import resample, Interpolation
    from isistools.csm2map.writers import write_geotiff

    model, body = load_camera("input.cub")
    grid = grid_from_params(crs=..., resolution=6.0, ...)
    coord_map = compute_transform_coarse(model, grid, body.radius_mean_m)
    projected = resample(data, coord_map)
    write_geotiff("output.tif", projected, grid)

Public API
----------
csm2map : function
    The main pipeline entry point.  Same name as the CLI command.
TargetBody : dataclass
    Target body ellipsoid (radii, NAIF ID), parsed from ALE's ISD.
OutputGrid : dataclass
    Output raster grid definition (CRS, affine, dimensions, bounds).
CoordinateMap : dataclass
    Dense mapping from output-pixel to input-pixel coordinates.
Interpolation : enum
    Resampling method (nearest, bilinear, bicubic).
"""

from isistools.csm2map.camera import TargetBody
from isistools.csm2map.grid import OutputGrid
from isistools.csm2map.pipeline import csm2map
from isistools.csm2map.resample import Interpolation
from isistools.csm2map.tiled import csm2map_tiled
from isistools.csm2map.transform import CoordinateMap

__all__ = [
    "CoordinateMap",
    "Interpolation",
    "OutputGrid",
    "TargetBody",
    "csm2map",
    "csm2map_tiled",
]
