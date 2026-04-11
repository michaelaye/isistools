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

    # ZSTD compression is dramatically faster than LZW at similar or
    # better compression ratios for float32 image data. ZSTD levels
    # 1-3 complete in well under half the wall time of LZW on a
    # 100 MB float32 cube while producing slightly smaller output.
    # The NUM_THREADS=ALL_CPUS option lets GDAL parallelize the
    # compression across cores.
    profile = {
        "driver": "GTiff",
        "dtype": data.dtype,
        "width": width,
        "height": height,
        "count": n_bands,
        "crs": grid.crs.to_wkt(),
        "transform": grid.transform,
        "nodata": nodata,
        "compress": "zstd",
        "zstd_level": 3,
        "num_threads": "ALL_CPUS",
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
    }

    with rasterio.open(str(output_path), "w", **profile) as dst:
        for b in range(n_bands):
            dst.write(data[b], b + 1)

    return output_path
