"""CTX EDR ingestion — Python replacement for ISIS mroctx2isis.

Reads a PDS3 CTX EDR file (.IMG), applies SQROOT decompression,
extracts dark/buffer pixels, and returns the image as a numpy array
with associated metadata.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pvl

# SQROOT lookup table: maps 8-bit compressed DN -> 12-bit uncompressed DN.
# Built from $ISISDATA/mro/calibration/ctxsqroot_001.lut
# Format: index = compressed value, value = decompressed value
_SQROOT_LUT = np.array(
    [
        1,
        3,
        5,
        7,
        9,
        11,
        13,
        15,
        17,
        20,
        22,
        24,
        27,
        29,
        32,
        35,
        38,
        41,
        44,
        47,
        50,
        54,
        58,
        61,
        65,
        69,
        73,
        77,
        82,
        86,
        91,
        95,
        100,
        105,
        110,
        115,
        121,
        126,
        131,
        137,
        143,
        149,
        155,
        161,
        167,
        173,
        179,
        186,
        193,
        199,
        206,
        213,
        220,
        228,
        235,
        243,
        250,
        258,
        266,
        274,
        282,
        290,
        298,
        306,
        315,
        324,
        332,
        341,
        350,
        359,
        369,
        378,
        387,
        397,
        407,
        416,
        426,
        436,
        446,
        457,
        467,
        478,
        488,
        499,
        510,
        521,
        532,
        543,
        554,
        566,
        577,
        589,
        601,
        613,
        625,
        637,
        649,
        662,
        674,
        687,
        699,
        712,
        725,
        738,
        751,
        765,
        778,
        792,
        805,
        819,
        833,
        847,
        861,
        875,
        890,
        904,
        919,
        933,
        948,
        963,
        978,
        993,
        1009,
        1024,
        1039,
        1055,
        1071,
        1087,
        1103,
        1119,
        1135,
        1151,
        1168,
        1184,
        1201,
        1218,
        1235,
        1252,
        1269,
        1286,
        1304,
        1321,
        1339,
        1356,
        1374,
        1392,
        1410,
        1429,
        1447,
        1465,
        1484,
        1502,
        1521,
        1540,
        1559,
        1578,
        1598,
        1617,
        1636,
        1656,
        1676,
        1696,
        1715,
        1736,
        1756,
        1776,
        1796,
        1817,
        1838,
        1858,
        1879,
        1900,
        1921,
        1943,
        1964,
        1985,
        2007,
        2029,
        2050,
        2072,
        2094,
        2117,
        2139,
        2161,
        2184,
        2206,
        2229,
        2252,
        2275,
        2298,
        2321,
        2345,
        2368,
        2392,
        2415,
        2439,
        2463,
        2487,
        2511,
        2535,
        2560,
        2584,
        2609,
        2634,
        2658,
        2683,
        2709,
        2734,
        2759,
        2784,
        2810,
        2836,
        2861,
        2887,
        2913,
        2939,
        2966,
        2992,
        3019,
        3045,
        3072,
        3099,
        3126,
        3153,
        3180,
        3207,
        3235,
        3262,
        3290,
        3317,
        3345,
        3373,
        3401,
        3430,
        3458,
        3486,
        3515,
        3544,
        3573,
        3601,
        3630,
        3660,
        3689,
        3718,
        3748,
        3777,
        3807,
        3837,
        3867,
        3897,
        3927,
        3958,
        3988,
        4019,
        4049,
        4080,
    ],
    dtype=np.int16,
)


@dataclass
class CTXMetadata:
    """Metadata extracted from a CTX EDR PDS3 label."""

    product_id: str
    spacecraft_clock_count: str
    start_time: str
    target_name: str
    mission_phase_name: str

    # Detector configuration
    spatial_summing: int
    sample_first_pixel: int
    line_exposure_duration: float  # milliseconds
    focal_plane_temperature: float  # Kelvin
    offset_mode_id: str
    sample_bit_mode_id: str

    # Image dimensions (from PDS3 IMAGE object)
    lines: int
    line_samples: int  # total samples per line including prefix/dark

    # Archive
    data_set_id: str
    orbit_number: int

    # Dark pixel configuration (computed from summing/edit mode)
    dark_start: int = 0  # first dark pixel index in raw line
    dark_end: int = 0  # last dark pixel index in raw line (inclusive)
    n_image_samples: int = 0  # number of science image samples

    # Raw dark pixels extracted during ingestion
    dark_pixels: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=np.int16))


def _parse_metadata(label: pvl.PVLModule) -> CTXMetadata:
    """Extract CTX metadata from a PDS3 label."""
    image = label["IMAGE"]

    # Spatial summing can be either keyword
    spatial_summing = int(label.get("SPATIAL_SUMMING", label.get("SAMPLING_FACTOR", 1)))

    # Sample first pixel (edit mode) can be either keyword
    sample_first_pixel = int(label.get("SAMPLE_FIRST_PIXEL", label.get("EDIT_MODE_ID", 0)))

    meta = CTXMetadata(
        product_id=str(label["PRODUCT_ID"]),
        spacecraft_clock_count=str(
            label.get(
                "SPACECRAFT_CLOCK_START_COUNT",
                label.get("SPACECRAFT_CLOCK_COUNT", ""),
            )
        ),
        start_time=str(label.get("START_TIME", "")),
        target_name=str(label.get("TARGET_NAME", "Mars")),
        mission_phase_name=str(label.get("MISSION_PHASE_NAME", "")),
        spatial_summing=spatial_summing,
        sample_first_pixel=sample_first_pixel,
        line_exposure_duration=float(label.get("LINE_EXPOSURE_DURATION", 0.0)),
        focal_plane_temperature=float(label.get("FOCAL_PLANE_TEMPERATURE", 0.0)),
        offset_mode_id=str(label.get("OFFSET_MODE_ID", "")),
        sample_bit_mode_id=str(label.get("SAMPLE_BIT_MODE_ID", "SQROOT")),
        lines=int(image["LINES"]),
        line_samples=int(image["LINE_SAMPLES"]),
        data_set_id=str(label.get("DATA_SET_ID", "")),
        orbit_number=int(label.get("ORBIT_NUMBER", 0)),
    )

    # Compute dark pixel ranges based on summing mode and edit mode.
    # Follows ISIS mroctx2isis logic exactly.
    if spatial_summing == 1:
        if sample_first_pixel == 0:
            # Edit Mode 0, summing 1: prefix columns 14-37
            meta.dark_start = 14
            meta.dark_end = 37
        else:
            # Edit Mode > 0, summing 1: prefix columns 0-15
            meta.dark_start = 0
            meta.dark_end = 15
    else:
        # Summing mode 2
        if sample_first_pixel == 0:
            meta.dark_start = 7
            meta.dark_end = 18
        else:
            meta.dark_start = 0
            meta.dark_end = 7

    # Suffix pixels
    if sample_first_pixel == 0:
        suffix = 18 if spatial_summing == 1 else 9
    else:
        suffix = 0

    # Image samples = total line samples - dark end - suffix - 1
    meta.n_image_samples = meta.line_samples - meta.dark_end - suffix - 1

    return meta


def ingest_ctx_edr(
    edr_path: str | Path,
    *,
    fill_gap: bool = True,
    preserve_special: bool = False,
) -> tuple[np.ndarray, CTXMetadata]:
    """Read a CTX PDS3 EDR and return the decompressed image + metadata.

    Parameters
    ----------
    edr_path : path-like
        Path to the PDS3 .IMG file.
    fill_gap : bool
        If True (default), treat raw DN 0 as NULL.
        Matches ISIS mroctx2isis FILLGAP=true behavior.
    preserve_special : bool
        If True, use ISIS float32 sentinel values (NULL, LIS, HIS, etc.)
        instead of NaN for special pixels.  This preserves the distinction
        between data gaps (NULL), detector saturation (HIS/LIS), and
        representation limits (HRS/LRS) that ISIS downstream tools rely on.
        If False (default), all special pixels become NaN.

    Returns
    -------
    image : np.ndarray
        2D float32 array of shape (lines, n_image_samples) with SQROOT-
        decompressed DN values.
    metadata : CTXMetadata
        Extracted metadata including dark pixel array.
    """
    from isistools.special_pixels import HIS, NULL

    edr_path = Path(edr_path)
    label = pvl.load(str(edr_path))
    meta = _parse_metadata(label)

    # Validate
    if meta.data_set_id and "EDR" not in meta.data_set_id:
        raise ValueError(f"Not a CTX EDR: DATA_SET_ID = {meta.data_set_id}")
    if meta.sample_bit_mode_id != "SQROOT":
        raise ValueError(f"Only SQROOT encoding is supported, got {meta.sample_bit_mode_id}")

    # Read raw image data.
    # PDS3 ^IMAGE pointer gives the starting record (1-based).
    image_pointer = label.get("^IMAGE", 2)
    if isinstance(image_pointer, pvl.Units):
        image_pointer = image_pointer.value
    record_bytes = int(label["RECORD_BYTES"])
    start_byte = (int(image_pointer) - 1) * record_bytes

    # Read all lines as uint8
    raw_bytes = np.fromfile(
        str(edr_path),
        dtype=np.uint8,
        offset=start_byte,
        count=meta.lines * meta.line_samples,
    ).reshape(meta.lines, meta.line_samples)

    # Extract dark pixels before SQROOT decompression
    dark_raw = raw_bytes[:, meta.dark_start : meta.dark_end + 1].copy()

    # Apply SQROOT decompression via LUT
    dark_decompressed = _SQROOT_LUT[dark_raw].astype(np.int16)

    # Handle fill_gap for dark pixels
    if fill_gap:
        dark_decompressed[dark_raw == 0] = -32768  # ISIS NULL for int16

    meta.dark_pixels = dark_decompressed

    # Extract image pixels (everything after the dark/prefix region)
    image_start = meta.dark_end + 1
    image_end = image_start + meta.n_image_samples
    image_raw = raw_bytes[:, image_start:image_end]

    # Apply SQROOT decompression
    image = _SQROOT_LUT[image_raw].astype(np.float32)

    # Map 8-bit special pixel values to float32 sentinels or NaN.
    # ISIS mroctx2isis maps: 0 → NULL (gap), 255 → HIS (high saturation).
    # LIS (low saturation) is also 0 in the 8-bit domain but only
    # distinguishable from NULL by context — we follow ISIS and treat
    # raw 0 as NULL when fill_gap is set.
    null_val = NULL if preserve_special else np.float32(np.nan)
    his_val = HIS if preserve_special else np.float32(np.nan)

    if fill_gap:
        image[image_raw == 0] = null_val
    # Raw 255 = high instrument saturation (max 8-bit value)
    image[image_raw == 255] = his_val

    return image, meta
