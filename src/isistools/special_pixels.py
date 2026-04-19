"""ISIS-compatible special pixel handling.

Provides constants and utilities for the 5 ISIS special pixel types,
both as float32 (for working arrays) and as sentinel values that can
be tested with standard numpy operations.

The special pixel types carry distinct meanings:

- **NULL**: No data — pixel was never acquired (gap in telemetry)
- **LIS**: Low Instrument Saturation — detector saturated low
- **LRS**: Low Representation Saturation — value below valid range
- **HIS**: High Instrument Saturation — detector saturated high
- **HRS**: High Representation Saturation — value above valid range

By default, ctxpipe/hirisepipe map all special pixels to NaN for
easy processing with numpy/scipy (which understand NaN natively).
Set ``preserve_special=True`` to keep the ISIS float32 sentinel values
instead, for workflows that need to distinguish between pixel types.
"""

from __future__ import annotations

import numpy as np

# IEEE 754 float32 sentinel values matching ISIS3 SpecialPixel.h.
# These are large negative numbers that are distinguishable from each
# other and from any valid pixel value.  They are ordered so that
# NULL < LRS < LIS < HRS < HIS < VALID_MIN.
NULL = np.float32(-3.4028226550889045e38)
LRS = np.float32(-3.4028228579130005e38)
LIS = np.float32(-3.4028230607370965e38)
HRS = np.float32(-3.4028232635611926e38)
HIS = np.float32(-3.4028234663852886e38)
# VALID_MIN is the boundary: anything <= this is a special pixel.
# Use float64 for the comparison threshold to avoid float32 overflow.
VALID_MIN = np.float64(-3.4028235677973366e38)

# 8-bit (unsigned char) special pixel values used in raw PDS data
NULL_U8 = np.uint8(0)
LIS_U8 = np.uint8(0)  # some instruments use 0 for LIS
HIS_U8 = np.uint8(255)  # some instruments use 255 for HIS

# 16-bit (SignedWord) special pixel values used in ISIS cubes
NULL_I16 = np.int16(-32768)
LRS_I16 = np.int16(-32767)
LIS_I16 = np.int16(-32766)
HIS_I16 = np.int16(-32765)
HRS_I16 = np.int16(-32764)
VALID_MIN_I16 = np.int16(-32752)


def is_special(pixel: np.ndarray) -> np.ndarray:
    """Test whether pixels are ISIS special values (any type).

    Works for both float32 (ISIS sentinel values) and int16 (ISIS cube
    SignedWord values).

    Parameters
    ----------
    pixel : np.ndarray
        Array of pixel values.

    Returns
    -------
    np.ndarray of bool
        True where pixel is any special type.
    """
    if pixel.dtype in (np.float32, np.float64):
        return ~np.isfinite(pixel) | (pixel <= float(VALID_MIN))
    elif pixel.dtype in (np.int16,):
        return pixel <= VALID_MIN_I16
    else:
        return np.zeros(pixel.shape, dtype=bool)


def special_to_nan(image: np.ndarray) -> np.ndarray:
    """Replace all ISIS special pixels with NaN.

    Parameters
    ----------
    image : np.ndarray
        Image with ISIS special pixel sentinels.

    Returns
    -------
    np.ndarray
        Float32 image with NaN replacing all special pixels.
    """
    out = image.astype(np.float32, copy=True)
    out[is_special(image)] = np.nan
    return out


def nan_to_special(image: np.ndarray, special_type: np.float32 = NULL) -> np.ndarray:
    """Replace NaN with a specific ISIS special pixel value.

    Parameters
    ----------
    image : np.ndarray
        Image with NaN for missing data.
    special_type : float32
        Which special pixel to use (default NULL).

    Returns
    -------
    np.ndarray
        Image with ISIS sentinel values instead of NaN.
    """
    out = image.astype(np.float32, copy=True)
    out[np.isnan(out)] = special_type
    return out
