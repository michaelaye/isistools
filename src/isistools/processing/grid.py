"""Define the output map-projected grid.

Supports two modes:
1. From an ISIS MAP PVL file (for cam2map compatibility)
2. From explicit projection + resolution + bounds parameters

The grid is defined as a rasterio-compatible affine transform + CRS + shape.
"""

import os
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pvl
from pyproj import CRS, Transformer

from isistools.geo.projections import mapping_to_crs

# pyproj raises an error when transforming between Earth (EPSG:4326) and
# non-Earth CRS (e.g. Mars projections). This is expected for planetary work.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")


@dataclass
class OutputGrid:
    """Fully specified output raster grid."""

    crs: CRS
    transform: "rasterio.transform.Affine"  # noqa: F821
    width: int  # n_samples
    height: int  # n_lines

    # Ground range in the projection's native units (usually meters)
    x_min: float
    x_max: float
    y_min: float
    y_max: float

    # Pixel resolution in meters/pixel
    resolution: float

    # Ground range in lat/lon (planetocentric, positive east)
    lat_min: float
    lat_max: float
    lon_min: float
    lon_max: float

    def pixel_to_map(self) -> tuple[np.ndarray, np.ndarray]:
        """Return 2D arrays of map (x, y) coordinates for all pixel centers."""
        cols = np.arange(self.width)
        rows = np.arange(self.height)
        cc, rr = np.meshgrid(cols, rows)
        # Apply affine: pixel center = transform * (col + 0.5, row + 0.5)
        x = self.transform.c + (cc + 0.5) * self.transform.a
        y = self.transform.f + (rr + 0.5) * self.transform.e
        return x, y

    def pixel_to_latlon(self) -> tuple[np.ndarray, np.ndarray]:
        """Return 2D arrays of (lat_rad, lon_rad) for all pixel centers.

        Uses pyproj to invert the map projection.
        Returns planetocentric lat/lon in radians.
        """
        x, y = self.pixel_to_map()
        # Inverse projection: map coords -> geographic lat/lon in degrees
        transformer = Transformer.from_crs(self.crs, CRS.from_epsg(4326), always_xy=True)
        lon_deg, lat_deg = transformer.transform(x, y)
        return np.deg2rad(lat_deg), np.deg2rad(lon_deg)


def grid_from_map_file(
    map_path: str | Path,
    camera_lat_range: tuple[float, float] | None = None,
    camera_lon_range: tuple[float, float] | None = None,
    resolution_override: float | None = None,
) -> OutputGrid:
    """Build an OutputGrid from an ISIS MAP PVL file.

    Parameters
    ----------
    map_path : path-like
        ISIS-style MAP file (.map) containing a Mapping group.
    camera_lat_range : tuple, optional
        (min_lat, max_lat) in degrees, from camera model.  Used if MAP file
        lacks range keywords.
    camera_lon_range : tuple, optional
        (min_lon, max_lon) in degrees, from camera model.
    resolution_override : float, optional
        Override pixel resolution in meters/pixel.

    Returns
    -------
    OutputGrid
    """
    label = pvl.load(str(map_path))
    mapping = label["Mapping"] if "Mapping" in label else label.get("Group", label)

    crs = mapping_to_crs(mapping)

    # Resolution
    if resolution_override is not None:
        res = resolution_override
    elif "PixelResolution" in mapping:
        res = float(mapping["PixelResolution"])  # meters/pixel
    elif "Scale" in mapping:
        # Scale is pixels/degree; convert via equatorial radius
        eq_r = float(mapping.get("EquatorialRadius", 3396190.0))
        scale = float(mapping["Scale"])
        res = (np.pi * eq_r) / (180.0 * scale)
    else:
        msg = "MAP file must contain PixelResolution or Scale"
        raise ValueError(msg)

    # If the MAP file includes an explicit UpperLeftCorner with Samples/Lines
    # (e.g. from an existing ISIS output cube), honor it verbatim. This
    # guarantees pixel-perfect grid alignment when replicating an ISIS run.
    has_upper_left = all(
        k in mapping for k in ("UpperLeftCornerX", "UpperLeftCornerY")
    )
    core = label.get("Core") if isinstance(label, dict) or hasattr(label, "get") else None
    n_samples = None
    n_lines = None
    if core is not None:
        dims = core.get("Dimensions") if hasattr(core, "get") else None
        if dims is not None:
            n_samples = int(dims.get("Samples", 0)) or None
            n_lines = int(dims.get("Lines", 0)) or None

    if has_upper_left and n_samples and n_lines:
        x_ul = float(mapping["UpperLeftCornerX"])
        y_ul = float(mapping["UpperLeftCornerY"])
        # Ground range (keep original lat/lon for reference)
        lat_min = float(mapping.get("MinimumLatitude", 0.0))
        lat_max = float(mapping.get("MaximumLatitude", 0.0))
        lon_min = float(mapping.get("MinimumLongitude", 0.0))
        lon_max = float(mapping.get("MaximumLongitude", 0.0))

        import rasterio.transform
        transform = rasterio.transform.Affine(
            res, 0.0, x_ul,
            0.0, -res, y_ul,
        )
        return OutputGrid(
            crs=crs,
            transform=transform,
            width=n_samples,
            height=n_lines,
            x_min=x_ul,
            x_max=x_ul + n_samples * res,
            y_min=y_ul - n_lines * res,
            y_max=y_ul,
            resolution=res,
            lat_min=lat_min,
            lat_max=lat_max,
            lon_min=lon_min,
            lon_max=lon_max,
        )

    # Ground range
    lat_min = _get_range(mapping, "MinimumLatitude", camera_lat_range, 0)
    lat_max = _get_range(mapping, "MaximumLatitude", camera_lat_range, 1)
    lon_min = _get_range(mapping, "MinimumLongitude", camera_lon_range, 0)
    lon_max = _get_range(mapping, "MaximumLongitude", camera_lon_range, 1)

    return grid_from_params(
        crs=crs,
        resolution=res,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
    )


def grid_from_params(
    crs: CRS | str,
    resolution: float,
    lat_min: float,
    lat_max: float,
    lon_min: float,
    lon_max: float,
) -> OutputGrid:
    """Build an OutputGrid from explicit parameters.

    Parameters
    ----------
    crs : pyproj CRS or proj string
        Output map projection.
    resolution : float
        Pixel resolution in meters/pixel.
    lat_min, lat_max, lon_min, lon_max : float
        Ground range in degrees (planetocentric, positive east).

    Returns
    -------
    OutputGrid
    """
    import rasterio.transform

    if isinstance(crs, str):
        crs = CRS.from_user_input(crs)

    # Project corner lat/lons to map coordinates
    transformer = Transformer.from_crs(CRS.from_epsg(4326), crs, always_xy=True)
    x_corners, y_corners = transformer.transform(
        [lon_min, lon_max, lon_min, lon_max],
        [lat_min, lat_min, lat_max, lat_max],
    )
    x_min = min(x_corners)
    x_max = max(x_corners)
    y_min = min(y_corners)
    y_max = max(y_corners)

    # Apply ISIS-compatible snap rule: floor() for minX, ceil() for maxY.
    # See ISIS ProjectionFactory::CreateForCube. This anchors the grid to
    # integer multiples of the pixel resolution from the projection origin,
    # matching what ISIS cam2map produces.
    eps = 1.0e-6
    if abs(x_min % resolution) > eps and abs(resolution - abs(x_min % resolution)) > eps:
        x_min = np.floor(x_min / resolution) * resolution
    if x_max < x_min + resolution:
        x_max = x_min + resolution
    if abs(y_max % resolution) > eps and abs(resolution - abs(y_max % resolution)) > eps:
        y_max = np.ceil(y_max / resolution) * resolution
    if y_min > y_max - resolution:
        y_min = y_max - resolution

    # ISIS rounds to nearest integer pixel count
    width = int((x_max - x_min) / resolution + 0.5)
    height = int((y_max - y_min) / resolution + 0.5)

    # Affine: pixel (0,0) has its upper-left corner at (x_min, y_max).
    # Build transform explicitly to avoid rasterio recomputing pixel size.
    transform = rasterio.transform.Affine(
        resolution, 0.0, x_min,
        0.0, -resolution, y_max,
    )

    return OutputGrid(
        crs=crs,
        transform=transform,
        width=width,
        height=height,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        resolution=resolution,
        lat_min=lat_min,
        lat_max=lat_max,
        lon_min=lon_min,
        lon_max=lon_max,
    )


def _get_range(mapping: dict, key: str, fallback: tuple | None, idx: int) -> float:
    """Get a range value from MAP file, falling back to camera range."""
    if key in mapping:
        return float(mapping[key])
    if fallback is not None:
        return fallback[idx]
    msg = f"MAP file missing {key} and no camera range provided"
    raise ValueError(msg)
