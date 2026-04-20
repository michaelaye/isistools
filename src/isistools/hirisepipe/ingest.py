"""HiRISE PDS EDR ingestion — Python replacement for ISIS hi2isis.

Reads a HiRISE PDS3 EDR file directly, applies LUT decompression,
and extracts image data + ancillary data (buffer pixels, dark pixels,
calibration image) needed for radiometric calibration.

No ISIS cube intermediate is needed — data goes straight from PDS EDR
to numpy arrays ready for hical.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pvl


@dataclass
class HiRISEEDR:
    """Data extracted from a HiRISE PDS3 EDR file."""

    # Science image (LUT-decompressed to 14-bit DN)
    image: np.ndarray  # (n_lines, n_samples), int16

    # Calibration image (LUT-decompressed)
    cal_image: np.ndarray  # (n_cal_lines, n_samples), int16

    # Buffer pixels (prefix) per science line
    buffer_pixels: np.ndarray  # (n_lines, 12), uint8→int16

    # Dark reference pixels (suffix) per science line
    dark_pixels: np.ndarray  # (n_lines, 16), uint8→int16

    # Calibration buffer/dark pixels
    cal_buffer: np.ndarray  # (n_cal_lines, 12)
    cal_dark: np.ndarray  # (n_cal_lines, 16)

    # Metadata from label
    product_id: str
    observation_id: str
    cpmm_number: int
    channel_number: int
    tdi: int
    binning: int
    scan_exposure_duration: float  # microseconds
    start_time: str
    spacecraft_clock_start_count: str
    target_name: str
    fpa_positive_y_temperature: float
    fpa_negative_y_temperature: float
    filter_name: str
    trim_lines: int
    n_lines: int
    n_samples: int
    lut_type: str


def _read_lut(label: pvl.PVLModule, file_path: Path) -> np.ndarray | None:
    """Read the 8-bit→14-bit inverse LUT from the EDR file.

    The LUT in the EDR maps 14-bit input → 8-bit output (16384 entries).
    We need the inverse: 8-bit → 14-bit (256 entries).
    """
    lut_info = label.get("LOOKUP_TABLE")
    if lut_info is None:
        return None

    n_rows = int(lut_info["ROWS"])  # 16384
    if n_rows != 16384:
        return None

    # Find the LUT location in the file
    # It follows the SCIENCE_CHANNEL_TABLE and CPMM_ENGINEERING_TABLE
    # We'll compute its offset from the label pointers
    # The LUT is referenced by ^LOOKUP_TABLE in the label
    lut_pointer = label.get("^LOOKUP_TABLE")
    if lut_pointer is None:
        return None

    if isinstance(lut_pointer, (list, tuple)):
        lut_offset = int(lut_pointer[0]) - 1  # byte offset, 0-based
    elif isinstance(lut_pointer, pvl.Quantity):
        lut_offset = int(lut_pointer.value) - 1
    else:
        lut_offset = int(lut_pointer) - 1

    # Read the forward LUT (14-bit → 8-bit, 16384 uint8 values)
    forward_lut = np.fromfile(str(file_path), dtype=np.uint8, offset=lut_offset, count=16384)

    # Build inverse LUT: for each 8-bit output value, find the 14-bit input
    # that maps to it. Use the midpoint of the range.
    inverse_lut = np.zeros(256, dtype=np.int16)
    for out_val in range(256):
        # Find all 14-bit inputs that map to this 8-bit output
        inputs = np.where(forward_lut == out_val)[0]
        if len(inputs) > 0:
            inverse_lut[out_val] = int(inputs[len(inputs) // 2])  # midpoint

    return inverse_lut


def _get_float(obj, key: str, default: float = 0.0) -> float:
    """Extract a float from a PVL object, handling Quantity types."""
    val = obj.get(key, default)
    if isinstance(val, pvl.Quantity):
        return float(val.value)
    return float(val)


def ingest_hirise_edr(edr_path: str | Path) -> HiRISEEDR:
    """Read a HiRISE PDS3 EDR and extract all data needed for calibration.

    Parameters
    ----------
    edr_path : path-like
        Path to the PDS3 EDR .IMG file.

    Returns
    -------
    HiRISEEDR
        Extracted image data, calibration data, and metadata.
    """
    edr_path = Path(edr_path)
    label = pvl.load(str(edr_path))

    # Extract metadata
    inst = label.get("INSTRUMENT_SETTING_PARAMETERS", {})
    time_params = label.get("TIME_PARAMETERS", {})
    temp_params = label.get("TEMPERATURE_PARAMETERS", {})
    image_info = label["IMAGE"]

    n_lines = int(image_info["LINES"])
    n_samples = int(image_info["LINE_SAMPLES"])
    sample_bits = int(image_info["SAMPLE_BITS"])
    prefix_bytes = int(image_info.get("LINE_PREFIX_BYTES", 18))
    suffix_bytes = int(image_info.get("LINE_SUFFIX_BYTES", 16))

    cal_info = label.get("CALIBRATION_IMAGE", {})
    n_cal_lines = int(cal_info.get("LINES", 0))

    cpmm = int(inst.get("MRO:CPMM_NUMBER", 0))
    channel = int(inst.get("MRO:CHANNEL_NUMBER", 0))
    tdi = int(inst.get("MRO:TDI", 128))
    binning = int(inst.get("MRO:BINNING", 1))
    trim_lines = int(inst.get("MRO:TRIM_LINES", 0))
    filter_name = str(inst.get("FILTER_NAME", "RED"))
    lut_type = str(inst.get("MRO:LOOKUP_TABLE_TYPE", "Stored"))

    start_time_raw = time_params.get("START_TIME", "")
    if hasattr(start_time_raw, "strftime"):
        start_time = start_time_raw.strftime("%Y-%m-%dT%H:%M:%S.%f")
    else:
        start_time = str(start_time_raw).replace("+00:00", "")

    sclk = str(time_params.get("SPACECRAFT_CLOCK_START_COUNT", ""))

    # Read LUT for decompression
    inverse_lut = _read_lut(label, edr_path)

    # Determine data locations from label pointers
    # The file layout is:
    #   Label → ScienChannelTable → LUT → CPMMEngTable →
    #   CalibrationLines (prefix+image+suffix) →
    #   ScienceLines (prefix+image+suffix) → GapTable

    cal_image_pointer = label.get("^CALIBRATION_IMAGE")
    image_pointer = label.get("^IMAGE")

    # Convert pointers to byte offsets
    def _pointer_to_offset(ptr):
        if isinstance(ptr, (list, tuple)):
            return int(ptr[0]) - 1
        elif isinstance(ptr, pvl.Quantity):
            return int(ptr.value) - 1
        elif ptr is not None:
            return int(ptr) - 1
        return None

    cal_offset = _pointer_to_offset(cal_image_pointer)
    img_offset = _pointer_to_offset(image_pointer)

    # Read calibration image data
    line_bytes = prefix_bytes + n_samples + suffix_bytes
    cal_image = np.empty((n_cal_lines, n_samples), dtype=np.int16)
    cal_buffer = np.empty((n_cal_lines, 12), dtype=np.int16)
    cal_dark = np.empty((n_cal_lines, 16), dtype=np.int16)

    with open(edr_path, "rb") as f:
        if cal_offset is not None:
            # Calibration lines have prefix (18 bytes) + image + suffix (16 bytes)
            # Prefix: 6 bytes line ID + 12 bytes buffer pixels
            for i in range(n_cal_lines):
                f.seek(cal_offset + i * line_bytes)
                line_data = np.frombuffer(f.read(line_bytes), dtype=np.uint8)

                # Buffer pixels: bytes 6-17 (12 pixels)
                cal_buffer[i] = line_data[6:18].astype(np.int16)

                # Image pixels: bytes 18 to 18+n_samples
                raw_pixels = line_data[prefix_bytes : prefix_bytes + n_samples]
                if inverse_lut is not None and sample_bits == 8:
                    cal_image[i] = inverse_lut[raw_pixels]
                else:
                    cal_image[i] = raw_pixels.astype(np.int16)

                # Dark reference pixels: last 16 bytes
                cal_dark[i] = line_data[prefix_bytes + n_samples :].astype(np.int16)

        # Read science image data
        image = np.empty((n_lines, n_samples), dtype=np.int16)
        buffer_pixels = np.empty((n_lines, 12), dtype=np.int16)
        dark_pixels = np.empty((n_lines, 16), dtype=np.int16)

        for i in range(n_lines):
            f.seek(img_offset + i * line_bytes)
            line_data = np.frombuffer(f.read(line_bytes), dtype=np.uint8)

            # Buffer pixels
            buffer_pixels[i] = line_data[6:18].astype(np.int16)

            # Image pixels
            raw_pixels = line_data[prefix_bytes : prefix_bytes + n_samples]
            if inverse_lut is not None and sample_bits == 8:
                image[i] = inverse_lut[raw_pixels]
            else:
                image[i] = raw_pixels.astype(np.int16)

            # Dark reference pixels
            dark_pixels[i] = line_data[prefix_bytes + n_samples :].astype(np.int16)

    return HiRISEEDR(
        image=image,
        cal_image=cal_image,
        buffer_pixels=buffer_pixels,
        dark_pixels=dark_pixels,
        cal_buffer=cal_buffer,
        cal_dark=cal_dark,
        product_id=str(label.get("PRODUCT_ID", "")),
        observation_id=str(label.get("OBSERVATION_ID", "")),
        cpmm_number=cpmm,
        channel_number=channel,
        tdi=tdi,
        binning=binning,
        scan_exposure_duration=_get_float(inst, "MRO:SCAN_EXPOSURE_DURATION"),
        start_time=start_time,
        spacecraft_clock_start_count=sclk,
        target_name=str(label.get("TARGET_NAME", "Mars")),
        fpa_positive_y_temperature=_get_float(temp_params, "MRO:FPA_POSITIVE_Y_TEMPERATURE"),
        fpa_negative_y_temperature=_get_float(temp_params, "MRO:FPA_NEGATIVE_Y_TEMPERATURE"),
        filter_name=filter_name,
        trim_lines=trim_lines,
        n_lines=n_lines,
        n_samples=n_samples,
        lut_type=lut_type,
    )
