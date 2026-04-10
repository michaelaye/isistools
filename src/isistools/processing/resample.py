"""Resample input image using a precomputed coordinate map.

Uses scipy.ndimage.map_coordinates for the actual interpolation.
This is the component that would benefit most from GPU acceleration later
(e.g., cupy.ndimage.map_coordinates or cv2.remap with CUDA).
"""

from enum import Enum

import numpy as np
from scipy.ndimage import map_coordinates

from isistools.processing.transform import CoordinateMap


class Interpolation(str, Enum):
    """Interpolation method for resampling."""

    NEAREST = "nearest"
    BILINEAR = "bilinear"
    BICUBIC = "bicubic"


# Map our names to scipy order parameter
_INTERP_ORDER = {
    Interpolation.NEAREST: 0,
    Interpolation.BILINEAR: 1,
    Interpolation.BICUBIC: 3,
}


def resample(
    input_data: np.ndarray,
    coord_map: CoordinateMap,
    interpolation: Interpolation = Interpolation.BICUBIC,
    fill_value: float = 0.0,
) -> np.ndarray:
    """Resample input image data using a coordinate map.

    Parameters
    ----------
    input_data : ndarray, shape (n_lines, n_samples) or (n_bands, n_lines, n_samples)
        Input image pixel data.
    coord_map : CoordinateMap
        Precomputed mapping from output pixels to input pixels.
    interpolation : Interpolation
        Interpolation method.
    fill_value : float
        Value for output pixels that map outside the input image.

    Returns
    -------
    ndarray
        Resampled image with shape (height, width) or (n_bands, height, width).
    """
    order = _INTERP_ORDER[interpolation]

    if input_data.ndim == 2:
        return _resample_band(input_data, coord_map, order, fill_value)

    # Multi-band: resample each band
    n_bands = input_data.shape[0]
    h, w = coord_map.shape
    output = np.full((n_bands, h, w), fill_value, dtype=np.float32)

    for b in range(n_bands):
        output[b] = _resample_band(input_data[b], coord_map, order, fill_value)

    return output


def _resample_band(
    band_data: np.ndarray,
    coord_map: CoordinateMap,
    order: int,
    fill_value: float,
) -> np.ndarray:
    """Resample a single band."""
    # map_coordinates expects coordinates as (row, col) = (line, sample)
    coordinates = np.array(
        [
            coord_map.input_lines,
            coord_map.input_samples,
        ]
    )

    result = map_coordinates(
        band_data.astype(np.float64),
        coordinates,
        order=order,
        mode="constant",
        cval=fill_value,
    ).astype(np.float32)

    # Apply validity mask
    result[~coord_map.valid] = fill_value

    return result
