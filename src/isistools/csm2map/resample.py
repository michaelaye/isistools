"""Resample input image using a precomputed coordinate map.

Uses ``scipy.ndimage.map_coordinates`` for the actual interpolation,
split across CPU threads (``scipy`` releases the GIL in its C
implementation, so threading gives close to linear speedup on the
resampling stage, which is usually the wall-time-dominant step for
large outputs).
"""

import os
from concurrent.futures import ThreadPoolExecutor
from enum import Enum

import numpy as np
from scipy.ndimage import map_coordinates

from isistools.csm2map.transform import CoordinateMap


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
    workers: int | None = None,
) -> np.ndarray:
    """Resample a single band.

    Implementation notes for performance:

    - ``scipy.ndimage.map_coordinates`` releases the GIL in its C
      implementation, so we split the output image into horizontal
      stripes and process each in a worker thread. On this machine
      (CPython 3.12, scipy 1.14) threading gives ~2.5x speedup with
      4 workers — the resample is memory-bandwidth bound at some
      point but there's meaningful headroom on CPU parallelism.
    - Everything runs in ``float32``. scipy internally computes in
      double precision no matter what, so the dtype choice mostly
      affects the input/output copies; keeping float32 halves the
      memory traffic on those.
    - The ``coord_map`` line/sample arrays are already float32 coming
      out of ``_bilinear_upsample_pair``; the input cube is float32
      coming out of ``read_isis_cube_raw``. Both casts are no-ops on
      the fast path.
    """
    h, w = coord_map.shape

    if band_data.dtype != np.float32:
        band_data = band_data.astype(np.float32, copy=False)

    if workers is None:
        # One thread per physical core is close to optimal; using every
        # logical core (HT) gives no extra win and starts hurting memory
        # bandwidth. Half the CPU count is a reasonable compromise.
        workers = max(1, (os.cpu_count() or 1) // 2)

    # Stripes need to be large enough that thread-setup overhead is
    # negligible. For outputs below ~1 M pixels a single-threaded path
    # is faster than any threading overhead.
    if workers <= 1 or h * w < 1_048_576:
        # coord_map.coords is already a (2, h, w) float32 buffer in
        # production (coarse path). For the rare float64 case
        # (validation, tests) scipy's C code accepts both — let it
        # decide rather than allocating a fresh stack.
        result = map_coordinates(
            band_data,
            coord_map.coords,
            order=order,
            mode="constant",
            cval=fill_value,
            output=np.float32,
        )
    else:
        # Split the output image into ``workers`` contiguous row
        # stripes. Each stripe gets its own slice of the coord map
        # and produces its own slice of the output. Threads share
        # the input ``band_data`` read-only, which is fine — scipy's
        # C code doesn't mutate it.
        result = np.empty((h, w), dtype=np.float32)
        stripe = (h + workers - 1) // workers

        def _process_stripe(i: int) -> None:
            r0 = i * stripe
            r1 = min(r0 + stripe, h)
            if r0 >= r1:
                return
            # View into the consolidated (2, h, w) buffer — no copy,
            # no allocation. Each stripe sees its own row band of both
            # channels.
            coordinates = coord_map.coords[:, r0:r1]
            map_coordinates(
                band_data,
                coordinates,
                order=order,
                mode="constant",
                cval=fill_value,
                output=result[r0:r1],
            )

        with ThreadPoolExecutor(max_workers=workers) as pool:
            list(pool.map(_process_stripe, range(workers)))

    # Apply validity mask. The fast path handles the common
    # ``fill_value == 0`` case with a single in-place multiply
    # (bool broadcast as 0/1), which has zero transient. The fallback
    # accepts a one-shot ~(h,w) bool transient for non-zero fills.
    if fill_value == 0:
        result *= coord_map.valid
    else:
        np.copyto(result, fill_value, where=~coord_map.valid)

    return result
