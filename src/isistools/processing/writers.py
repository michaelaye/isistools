"""Output writers for map-projected data."""

from pathlib import Path

import numpy as np
import rasterio

from isistools.processing.grid import OutputGrid


def write_geotiff(
    output_path: str | Path,
    data: np.ndarray,
    grid: OutputGrid,
    nodata: float = 0.0,
) -> Path:
    """Write map-projected data as a GeoTIFF.

    Parameters
    ----------
    output_path : path-like
        Output file path (should end in .tif).
    data : ndarray, shape (height, width) or (n_bands, height, width)
        Map-projected pixel data.
    grid : OutputGrid
        Grid definition with CRS and affine transform.
    nodata : float
        NoData value.

    Returns
    -------
    Path to written file.
    """
    output_path = Path(output_path)

    if data.ndim == 2:
        data = data[np.newaxis, ...]
    n_bands, height, width = data.shape

    assert height == grid.height, f"Data height {height} != grid height {grid.height}"
    assert width == grid.width, f"Data width {width} != grid width {grid.width}"

    # Replace NaN with nodata
    data = np.where(np.isnan(data), nodata, data)

    profile = {
        "driver": "GTiff",
        "dtype": data.dtype,
        "width": width,
        "height": height,
        "count": n_bands,
        "crs": grid.crs.to_wkt(),
        "transform": grid.transform,
        "nodata": nodata,
        "compress": "lzw",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(str(output_path), "w", **profile) as dst:
        for b in range(n_bands):
            dst.write(data[b], b + 1)

    return output_path
