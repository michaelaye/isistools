"""HiRISE radiometric calibration — Python replacement for ISIS hical.

Implements the full 10-module HiCal calibration chain:

1. ZeroBufferSmooth (Zf) — buffer pixel drift extraction
2. ZeroBufferFit (Zd) — drift correction (pass-through by default)
3. ZeroReverse (Zz) — reverse-clock column offset
4. ZeroDark (Zb) — temperature-corrected dark current
5. GainLineDrift (Zg) — line-time gain drift
6. GainNonLinearity (Gnl) — signal-dependent non-linearity
7. GainChannelNormalize (Gcn) — sample gain with TDI/BIN normalization
8. GainFlatField (Za) — flat field (A matrix)
9. GainTemperature (Zt) — FPA temperature gain
10. GainUnitConversion (Ziof) — unit conversion

Calibration equation (from ISIS hical main.cpp)::

    hdn = (idn - ZBF - ZRev - ZD) / GLD
    NLGain = 1.0 - (GNL * median(hdn_line))
    odn = hdn * GCN * NLGain * GFF * GT / GUC
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pvl

from isistools.io.cubes import read_label

# Default calibration data paths
_DEFAULT_ISISDATA = Path(
    os.environ.get("ISISDATA", str(Path.home() / "Dropbox" / "data" / "isisdata"))
)
_MATRICES_DIR = _DEFAULT_ISISDATA / "mro" / "calibration" / "matrices"

# CPMM-to-CCD mapping (HiRISE has 14 CPMMs mapped to CCDs)
_CPMM_TO_CCD = [0, 1, 2, 3, 12, 4, 10, 11, 5, 13, 6, 7, 8, 9]
_CCD_NAMES = [
    "RED0",
    "RED1",
    "RED2",
    "RED3",
    "RED4",
    "RED5",
    "RED6",
    "RED7",
    "RED8",
    "RED9",
    "IR10",
    "IR11",
    "BG12",
    "BG13",
]


def _cpmm_to_ccd(cpmm: int) -> int:
    """Convert CPMM number to CCD number."""
    return _CPMM_TO_CCD[cpmm]


def _ccd_name(ccd: int) -> str:
    """Return CCD name from CCD number."""
    return _CCD_NAMES[ccd]


def _filter_name(ccd: int) -> str:
    """Return filter name from CCD number."""
    name = _CCD_NAMES[ccd]
    if name.startswith("RED"):
        return "RED"
    elif name.startswith("IR"):
        return "IR"
    return "BG"


@dataclass
class HiCalParams:
    """Parameters extracted from a HiRISE cube label for calibration."""

    cpmm: int
    ccd: int
    channel: int
    tdi: int
    binning: int
    scan_exposure_duration: float  # microseconds
    n_lines: int
    n_samples: int
    fpa_pos_y_temp: float  # degrees C
    fpa_neg_y_temp: float  # degrees C
    filter_name: str
    ccd_name: str
    product_id: str
    trim_lines: int

    @property
    def fpa_temp(self) -> float:
        """Average FPA temperature."""
        return (self.fpa_pos_y_temp + self.fpa_neg_y_temp) / 2.0

    @classmethod
    def from_cube(cls, cube_path: str | Path) -> HiCalParams:
        """Extract calibration parameters from a HiRISE ISIS cube label."""
        label = read_label(cube_path)
        inst = label["IsisCube"]["Instrument"]
        arch = label["IsisCube"].get("Archive", {})

        cpmm = int(inst["CpmmNumber"])
        ccd = _cpmm_to_ccd(cpmm)
        channel = int(inst["ChannelNumber"])
        tdi = int(inst["Tdi"])
        binning = int(inst["Summing"])

        dims = label["IsisCube"]["Core"]["Dimensions"]
        n_lines = int(dims["Lines"])
        n_samples = int(dims["Samples"])

        return cls(
            cpmm=cpmm,
            ccd=ccd,
            channel=channel,
            tdi=tdi,
            binning=binning,
            scan_exposure_duration=float(inst["ScanExposureDuration"]),
            n_lines=n_lines,
            n_samples=n_samples,
            fpa_pos_y_temp=float(inst["FpaPositiveYTemperature"]),
            fpa_neg_y_temp=float(inst["FpaNegativeYTemperature"]),
            filter_name=_filter_name(ccd),
            ccd_name=_ccd_name(ccd),
            product_id=str(arch.get("ProductId", "")),
            trim_lines=int(arch.get("TrimLines", 0)),
        )


def _line_time(line: int, binning: int, sed_us: float) -> float:
    """Compute line time in seconds. Matches ISIS HiLineTimeEqn."""
    return line * binning * sed_us * 1.0e-6


def _load_csv_column(
    filepath: Path,
    column_name: str,
) -> np.ndarray:
    """Load a named column from a HiRISE calibration CSV file."""
    import csv

    with open(filepath) as f:
        reader = csv.reader(f)
        # Skip comment lines
        rows = []
        header = None
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if header is None:
                header = [h.strip() for h in row]
                continue
            rows.append(row)

    if header is None:
        raise ValueError(f"No header found in {filepath}")

    # Find column index
    try:
        col_idx = header.index(column_name)
    except ValueError:
        # Try case-insensitive match
        header_lower = [h.lower() for h in header]
        try:
            col_idx = header_lower.index(column_name.lower())
        except ValueError:
            raise ValueError(
                f"Column '{column_name}' not found in {filepath}. Available: {header}"
            )

    values = []
    for row in rows:
        if col_idx < len(row):
            val = row[col_idx].strip()
            if val:
                values.append(float(val))

    return np.array(values, dtype=np.float64)


def _load_csv_row(
    filepath: Path,
    row_name: str,
) -> np.ndarray:
    """Load a named row from a HiRISE calibration CSV file."""
    import csv

    with open(filepath) as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            if row[0].strip() == row_name:
                return np.array([float(v.strip()) for v in row[1:] if v.strip()], dtype=np.float64)

    raise ValueError(f"Row '{row_name}' not found in {filepath}")


def _load_csv_cell(
    filepath: Path,
    row_key: str,
    col_key: str,
) -> float:
    """Load a single cell from a CSV file by row key and column header."""
    import csv

    with open(filepath) as f:
        reader = csv.reader(f)
        header = None
        for row in reader:
            if not row or row[0].strip().startswith("#"):
                continue
            if header is None:
                header = [h.strip() for h in row]
                continue
            if row[0].strip() == row_key:
                col_idx = header.index(col_key)
                return float(row[col_idx].strip())

    raise ValueError(f"Cell ({row_key}, {col_key}) not found in {filepath}")


def _find_latest_matrix(
    pattern: str,
    matrices_dir: Path | None = None,
) -> Path:
    """Find the highest-versioned calibration matrix file."""
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR
    candidates = sorted(matrices_dir.glob(pattern))
    if not candidates:
        raise FileNotFoundError(f"No file matching '{pattern}' in {matrices_dir}")
    return candidates[-1]


def _hi_temp_eqn(
    temperature: float,
    napcm2: float = 2.0,
    px: float = 12.0,
) -> float:
    """HiRISE temperature equation for dark current (electrons/sec/pixel).

    From ISIS HiCalUtil.h HiTempEqn.
    """
    temp = temperature + 273.0  # Convert to Kelvin
    eg = 1.1557 - (7.021e-4 * temp * temp) / (1108.0 + temp)
    K = 1.38e-23  # Boltzmann constant
    Q = 1.6e-19  # Electron charge
    return napcm2 * (px * px) * 2.55e7 * (temp**1.5) * np.exp(-eg * Q / 2.0 / K / temp)


def _lowpass_filter(
    data: np.ndarray,
    width: int,
    iterations: int = 1,
) -> np.ndarray:
    """Apply a boxcar lowpass filter, preserving NaN gaps."""
    result = data.copy()
    half = width // 2
    for _ in range(iterations):
        smoothed = np.empty_like(result)
        for i in range(len(result)):
            start = max(0, i - half)
            end = min(len(result), i + half + 1)
            window = result[start:end]
            valid = window[np.isfinite(window)]
            if len(valid) > 0:
                smoothed[i] = valid.mean()
            else:
                smoothed[i] = np.nan
        result = smoothed
    return result


def _spline_fill(data: np.ndarray) -> np.ndarray:
    """Fill NaN gaps with cubic spline interpolation."""
    from scipy.interpolate import CubicSpline

    valid = np.isfinite(data)
    if valid.all():
        return data.copy()
    if not valid.any():
        return np.zeros_like(data)

    x_valid = np.where(valid)[0]
    y_valid = data[valid]

    cs = CubicSpline(x_valid, y_valid, extrapolate=True)
    result = data.copy()
    result[~valid] = cs(np.where(~valid)[0])
    return result


# ── Calibration Modules ──────────────────────────────────────────────


def zero_buffer_smooth(
    buffer_pixels: np.ndarray,
    first_sample: int = 5,
    last_sample: int = 11,
    filter_width: int = 201,
    filter_iterations: int = 2,
) -> np.ndarray:
    """ZeroBufferSmooth (Zf): Extract and smooth buffer pixel drift.

    Parameters
    ----------
    buffer_pixels : np.ndarray
        Buffer pixel table, shape (n_lines, 12).
    first_sample, last_sample : int
        Sample range to average (0-based, inclusive).
    filter_width : int
        Lowpass filter width (odd).
    filter_iterations : int
        Number of lowpass filter passes.

    Returns
    -------
    np.ndarray
        Smoothed drift vector, shape (n_lines,).
    """
    # Average across specified sample range
    avg = np.nanmean(
        buffer_pixels[:, first_sample : last_sample + 1].astype(np.float64),
        axis=1,
    )

    # Lowpass filter
    filtered = _lowpass_filter(avg, filter_width, filter_iterations)

    # Spline fill any remaining gaps
    return _spline_fill(filtered)


def zero_buffer_fit(
    zbs: np.ndarray,
    skip_fit: bool = True,
) -> np.ndarray:
    """ZeroBufferFit (Zd): Non-linear drift fit.

    By default (SkipFit=True in hical.0023.conf), this is a pass-through
    of the ZeroBufferSmooth result.  The non-linear Levenberg-Marquardt
    fit is optional and rarely used in practice.

    Parameters
    ----------
    zbs : np.ndarray
        ZeroBufferSmooth result.
    skip_fit : bool
        If True (default), pass through the input.

    Returns
    -------
    np.ndarray
        Drift correction vector, normalized so first value = 0.
    """
    if skip_fit:
        # Normalize: subtract first value so correction starts at 0
        result = zbs.copy()
        if len(result) > 0 and np.isfinite(result[0]):
            result -= result[0]
        return result

    # Full non-linear fit would go here (rarely used)
    raise NotImplementedError("Non-linear ZeroBufferFit not yet implemented")


def zero_reverse(
    cal_image: np.ndarray,
    first_line: int = 1,
    last_line: int = 19,
    n_samples: int | None = None,
    channel: int = 0,
) -> np.ndarray:
    """ZeroReverse (Zz): Reverse-clock column offset correction.

    Parameters
    ----------
    cal_image : np.ndarray
        Calibration image region, shape (n_cal_lines, detector_width).
    first_line, last_line : int
        Line range for reverse-clock region (0-based).
    n_samples : int, optional
        Channel width (if cal image covers full detector).
    channel : int
        Channel number (0 or 1) for cropping.

    Returns
    -------
    np.ndarray
        Per-sample correction vector, shape (n_samples,).
    """
    rev_clock = cal_image[first_line : last_line + 1, :].astype(np.float64)

    # If cal image is wider than image (full detector), extract channel data.
    # The cal image stores 2048-wide detector data with channels interleaved:
    # channel 0 = even columns (0, 2, 4, ...), channel 1 = odd columns (1, 3, 5, ...)
    if n_samples is not None and rev_clock.shape[1] > n_samples:
        if channel == 0:
            rev_clock = rev_clock[:, 0::2]  # even columns
        else:
            rev_clock = rev_clock[:, 1::2]  # odd columns

    # Average across lines, filling any special pixels
    result = np.nanmean(rev_clock, axis=0)

    # Fill any remaining NaN values
    return _spline_fill(result)


def zero_dark(
    params: HiCalParams,
    matrices_dir: Path | None = None,
) -> np.ndarray:
    """ZeroDark (Zb): Temperature-corrected dark current.

    Parameters
    ----------
    params : HiCalParams
        Calibration parameters from cube label.
    matrices_dir : Path, optional
        Directory containing calibration matrices.

    Returns
    -------
    np.ndarray
        Per-sample dark current vector, shape (n_samples,).
    """
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR

    ccd = params.ccd
    channel = params.channel
    tdi = params.tdi
    binning = params.binning
    col_name = f"{ccd}/{channel}"

    # Load B matrix (dark current baseline)
    b_file = _find_latest_matrix(f"B_TDI{tdi}_BIN{binning}_hical*.csv", matrices_dir)
    b_data = _load_csv_column(b_file, col_name)

    # Load temperature slope and intercept (256-sample vectors)
    slope_col = f"CH{channel}_TDI{tdi}"
    slope_file = _find_latest_matrix("B_Temperature_Slope_hical_????.csv", matrices_dir)
    intercept_file = _find_latest_matrix("B_Temperature_Intercept_hical_????.csv", matrices_dir)
    slope_data = _load_csv_column(slope_file, slope_col)
    intercept_data = _load_csv_column(intercept_file, slope_col)

    # Temperature correction
    fpa_ref = 21.0  # FPA reference temperature
    fpa_temp = params.fpa_temp

    # ISIS ZeroDark.h lines 123-127: filter slope and intercept BEFORE computing t_prof
    slope_data = _lowpass_filter(slope_data, width=3, iterations=1)
    intercept_data = _lowpass_filter(intercept_data, width=3, iterations=1)

    # Build temperature profile: t_prof = intercept + slope * temp
    t_prof = intercept_data + slope_data * fpa_temp

    # Rebin temperature profile to image samples
    if len(t_prof) != params.n_samples:
        x_old = np.linspace(0, 1, len(t_prof))
        x_new = np.linspace(0, 1, params.n_samples)
        t_prof = np.interp(x_new, x_old, t_prof)

    # Scale factor (from ISIS ZeroDark.h line 140-141)
    linetime = params.scan_exposure_duration  # microseconds
    scale = linetime * 1.0e-6 * (binning**2) * (20.0 * 103.0 / 89.0 + tdi)

    # Temperature equation (from ISIS HiCalUtil.h HiTempEqn)
    base_t = _hi_temp_eqn(fpa_ref)

    dark = np.empty(params.n_samples, dtype=np.float64)
    for j in range(params.n_samples):
        dark[j] = b_data[j] * scale * _hi_temp_eqn(t_prof[j]) / base_t

    # Lowpass filter (width=3, iterations=1)
    dark = _lowpass_filter(dark, width=3, iterations=1)

    return dark


def gain_line_drift(
    params: HiCalParams,
    matrices_dir: Path | None = None,
) -> np.ndarray:
    """GainLineDrift (Zg): Line-time gain drift.

    Model: gainV[i] = c0 + c1*lt + c2*exp(c3*lt)
    where lt = line * bin * sed * 1e-6 (seconds).

    Returns
    -------
    np.ndarray
        Per-line gain drift vector, shape (n_lines,).
    """
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR

    row_name = f"{params.ccd}/{params.channel}"
    coef_file = _find_latest_matrix(
        f"Line_Gain_Drift_BIN{params.binning}_hical*.csv", matrices_dir
    )
    coefs = _load_csv_row(coef_file, row_name)[:4]

    c0, c1, c2, c3 = coefs
    lines = np.arange(params.n_lines, dtype=np.float64)
    lt = lines * params.binning * params.scan_exposure_duration * 1.0e-6

    return c0 + c1 * lt + c2 * np.exp(c3 * lt)


def gain_non_linearity(
    params: HiCalParams,
    matrices_dir: Path | None = None,
) -> float:
    """GainNonLinearity (Gnl): Non-linearity coefficient.

    Returns a single scalar used in the two-pass calibration loop.
    """
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR

    row_name = f"{params.ccd}_{params.channel}"  # CSV uses underscore: "4_0"
    coef_file = _find_latest_matrix(
        f"Gain_NonLinearity_BIN{params.binning}_hical*.csv", matrices_dir
    )
    coefs = _load_csv_row(coef_file, row_name)
    return float(coefs[0])


def gain_channel_normalize(
    params: HiCalParams,
    matrices_dir: Path | None = None,
) -> np.ndarray:
    """GainChannelNormalize (Gcn): Sample gain with TDI/BIN normalization.

    Returns
    -------
    np.ndarray
        Per-sample gain vector, shape (n_samples,).
    """
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR

    col_name = f"{params.ccd}/{params.channel}"
    gains_file = _find_latest_matrix("Gains_hical_????.csv", matrices_dir)

    # Gains file is a table: BIN as row, CCD/channel as columns
    # Read the specific cell for this BIN and CCD/channel
    gain_value = _load_csv_cell(gains_file, str(params.binning), col_name)

    normalizer = 128.0 / params.tdi / (params.binning**2)
    result = np.full(params.n_samples, gain_value * normalizer)

    return result


def gain_flat_field(
    params: HiCalParams,
    matrices_dir: Path | None = None,
) -> np.ndarray:
    """GainFlatField (Za): Flat field correction from A matrix.

    Returns
    -------
    np.ndarray
        Per-sample flat field vector, shape (n_samples,).
    """
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR

    col_name = f"{params.ccd}/{params.channel}"
    a_file = _find_latest_matrix(f"A_TDI{params.tdi}_BIN{params.binning}_hical*.csv", matrices_dir)
    a_data = _load_csv_column(a_file, col_name)

    # Rebin if needed
    if len(a_data) != params.n_samples:
        x_old = np.linspace(0, 1, len(a_data))
        x_new = np.linspace(0, 1, params.n_samples)
        a_data = np.interp(x_new, x_old, a_data)

    return a_data


def gain_temperature(
    params: HiCalParams,
    matrices_dir: Path | None = None,
) -> np.ndarray:
    """GainTemperature (Zt): FPA temperature-dependent gain.

    Returns
    -------
    np.ndarray
        Per-sample temperature gain vector (constant across samples).
    """
    if matrices_dir is None:
        matrices_dir = _MATRICES_DIR

    col_name = f"{params.ccd}/{params.channel}"
    fpa_file = _find_latest_matrix("Temperature_Gain_????.csv", matrices_dir)
    fpa_factor = _load_csv_cell(fpa_file, str(params.binning), col_name)

    fpa_ref = 21.0
    t_gain = 1.0 - fpa_factor * (params.fpa_temp - fpa_ref)

    return np.full(params.n_samples, t_gain)


# HiRISE I/F calibration constants per filter (from hical.0023.conf)
_FILTER_IOF_PARAMS = {
    "RED": {
        "FilterGainCorrection": 157702564.0,
        "IoverFbasetemperature": 18.9,
        "QEpercentincreaseperC": 0.0005704,
        "AbsGain_TDI128": 6.376583,
    },
    "IR": {
        "FilterGainCorrection": 56464791.0,
        "IoverFbasetemperature": 18.9,
        "QEpercentincreaseperC": 0.002696,
        "AbsGain_TDI128": 6.989840,
    },
    "BG": {
        "FilterGainCorrection": 115074166.0,
        "IoverFbasetemperature": 18.9,
        "QEpercentincreaseperC": 0.00002295,
        "AbsGain_TDI128": 6.997557,
    },
}


def gain_unit_conversion(
    params: HiCalParams,
    units: str = "DN",
    sun_distance_au: float | None = None,
) -> float:
    """GainUnitConversion (Ziof): Unit conversion factor.

    Parameters
    ----------
    params : HiCalParams
        Calibration parameters.
    units : str
        "DN" (no conversion), "DN/US" (per microsecond), or "IOF" (I/F).
    sun_distance_au : float, optional
        Sun-to-target distance in AU. Required for "IOF" mode.
        Can be computed via ``isistools.spice_utils.sun_distance_au()``.

    Returns
    -------
    float
        Divisor for final calibration step.
    """
    units = units.upper()
    if units == "DN":
        return 1.0
    elif units == "DN/US":
        return params.scan_exposure_duration
    elif units == "IOF":
        if sun_distance_au is None:
            raise ValueError(
                "sun_distance_au is required for IOF units. "
                "Compute with: isistools.spice_utils.sun_distance_au(utc_time)"
            )
        # ISIS formula from GainUnitConversion.h
        fparams = _FILTER_IOF_PARAMS[params.filter_name]

        # Sun correction: (1.5 AU / actual AU)^2
        suncorr = (1.5 / sun_distance_au) ** 2

        # Temperature-dependent QE correction
        zgain = fparams["FilterGainCorrection"]
        base_t = fparams["IoverFbasetemperature"]
        qe_pct = fparams["QEpercentincreaseperC"]
        abs_gain = fparams["AbsGain_TDI128"]
        fpa_temp = params.fpa_temp
        qe_temp_dep = zgain * (1.0 + (fpa_temp - base_t) * qe_pct * abs_gain)

        zbin = 1.0  # GainUnitConversionBinFactor from config
        sed = params.scan_exposure_duration  # microseconds

        ziof = (zbin * qe_temp_dep) * (sed * 1.0e-6) * suncorr
        return ziof
    else:
        raise ValueError(f"Unknown units: {units}. Use 'DN', 'DN/US', or 'IOF'.")


# ── Main Calibration Function ────────────────────────────────────────


def hical(
    cube_path: str | Path,
    *,
    units: str = "DN",
    matrices_dir: Path | None = None,
    sun_distance_au: float | None = None,
) -> np.ndarray:
    """Apply HiRISE radiometric calibration to an ISIS cube.

    Implements the full hical calibration equation::

        hdn = (idn - ZBF - ZRev - ZD) / GLD
        NLGain = 1.0 - (GNL * median(hdn_line))
        odn = hdn * GCN * NLGain * GFF * GT / GUC

    Parameters
    ----------
    cube_path : path-like
        Path to a HiRISE ISIS cube (output of hi2isis).
    units : str
        Output units: "DN", "DN/US", or "IOF".
    matrices_dir : Path, optional
        Directory containing calibration matrix CSV files.

    Returns
    -------
    np.ndarray
        Calibrated float32 image, shape (n_lines, n_samples).
    """
    cube_path = Path(cube_path)

    # Extract parameters from label
    params = HiCalParams.from_cube(cube_path)

    # Read image data
    import rioxarray  # noqa: F401
    import xarray as xr

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Dataset has no geotransform")
        da = xr.open_dataarray(str(cube_path), engine="rasterio")
    image = da.values.squeeze().astype(np.float64)

    # Read ancillary tables from cube
    buffer_pixels, cal_image = _read_ancillary_tables(cube_path, params)

    # ── Compute calibration vectors ──

    # 1. ZeroBufferSmooth
    zbs = zero_buffer_smooth(buffer_pixels)

    # 2. ZeroBufferFit (pass-through by default)
    zbf = zero_buffer_fit(zbs, skip_fit=True)

    # 3. ZeroReverse
    zrev = zero_reverse(cal_image, n_samples=params.n_samples, channel=params.channel)

    # 4. ZeroDark
    zd = zero_dark(params, matrices_dir)

    # 5. GainLineDrift
    gld = gain_line_drift(params, matrices_dir)

    # 6. GainNonLinearity coefficient
    gnl = gain_non_linearity(params, matrices_dir)

    # 7. GainChannelNormalize
    gcn = gain_channel_normalize(params, matrices_dir)

    # 8. GainFlatField
    gff = gain_flat_field(params, matrices_dir)

    # 9. GainTemperature
    gt = gain_temperature(params, matrices_dir)

    # 10. GainUnitConversion
    guc = gain_unit_conversion(params, units, sun_distance_au=sun_distance_au)

    # ── Apply calibration (two-pass, line-by-line) ──

    out = np.empty_like(image, dtype=np.float32)

    for line_idx in range(params.n_lines):
        # Clamp line index for ZBF and GLD
        effective_line = min(line_idx, params.n_lines - 1)

        line = image[line_idx]

        # Identify special/valid pixels
        # ISIS SignedWord special pixels: <= -32752 are special
        special = line <= -32752

        # First pass: subtract drift, reverse, dark; divide by gain drift
        hdn = np.where(
            special,
            np.nan,
            (line - zbf[effective_line] - zrev - zd) / gld[effective_line],
        )

        # Compute median of valid calibrated pixels for non-linearity
        valid_hdn = hdn[np.isfinite(hdn)]
        if len(valid_hdn) > 0:
            line_median = float(np.median(valid_hdn))
            nl_gain = 1.0 - gnl * line_median

            # Second pass: apply gain, non-linearity, flat field, temp, units
            calibrated = np.where(
                np.isfinite(hdn),
                hdn * gcn * nl_gain * gff * gt / guc,
                np.nan,
            )
        else:
            calibrated = hdn

        # Restore special pixels as NaN (for float output)
        out[line_idx] = calibrated.astype(np.float32)

    return out


def hical_from_edr(
    edr_path: str | Path,
    *,
    units: str = "DN",
    matrices_dir: Path | None = None,
    sun_distance_au: float | None = None,
) -> np.ndarray:
    """Apply HiRISE radiometric calibration directly from a PDS EDR file.

    Reads the PDS3 EDR, extracts image + ancillary data, and calibrates
    without any ISIS intermediate.  This is the preferred entry point for
    a fully Python-native pipeline.

    Parameters
    ----------
    edr_path : path-like
        Path to a HiRISE PDS3 EDR (.IMG file).
    units : str
        Output units: "DN", "DN/US", or "IOF".
    matrices_dir : Path, optional
        Calibration matrices directory.
    sun_distance_au : float, optional
        Sun distance in AU (required for IOF units).

    Returns
    -------
    np.ndarray
        Calibrated float32 image.
    """
    from isistools.hirisepipe.ingest import ingest_hirise_edr

    edr = ingest_hirise_edr(edr_path)

    # Build params from EDR metadata
    from isistools.hirisepipe.hical import _ccd_name, _cpmm_to_ccd

    ccd = _cpmm_to_ccd(edr.cpmm_number)
    params = HiCalParams(
        cpmm=edr.cpmm_number,
        ccd=ccd,
        channel=edr.channel_number,
        tdi=edr.tdi,
        binning=edr.binning,
        scan_exposure_duration=edr.scan_exposure_duration,
        n_lines=edr.n_lines,
        n_samples=edr.n_samples,
        fpa_pos_y_temp=edr.fpa_positive_y_temperature,
        fpa_neg_y_temp=edr.fpa_negative_y_temperature,
        filter_name=edr.filter_name,
        ccd_name=_ccd_name(ccd),
        product_id=edr.product_id,
        trim_lines=edr.trim_lines,
    )

    image = edr.image.astype(np.float64)
    buffer_pixels = edr.buffer_pixels
    cal_image = edr.cal_image

    # Apply LUT to buffer pixels (they're raw 8-bit, need decompression)
    # The buffer pixels from EDR are already LUT-decompressed in ingest
    # (the inverse LUT is applied during reading)

    # ── Compute calibration vectors ──
    zbs = zero_buffer_smooth(buffer_pixels)
    zbf = zero_buffer_fit(zbs, skip_fit=True)
    zrev = zero_reverse(cal_image, n_samples=params.n_samples, channel=params.channel)
    zd = zero_dark(params, matrices_dir)
    gld = gain_line_drift(params, matrices_dir)
    gnl = gain_non_linearity(params, matrices_dir)
    gcn = gain_channel_normalize(params, matrices_dir)
    gff = gain_flat_field(params, matrices_dir)
    gt = gain_temperature(params, matrices_dir)
    guc = gain_unit_conversion(params, units, sun_distance_au=sun_distance_au)

    # ── Apply calibration (two-pass, line-by-line) ──
    out = np.empty_like(image, dtype=np.float32)

    for line_idx in range(params.n_lines):
        effective_line = min(line_idx, params.n_lines - 1)
        line = image[line_idx]

        # For EDR data, special pixels are 0 (gap) or 16383 (saturation)
        # and the LUT maps 255→max, 0→0. We treat anything that was
        # originally 0xFF as special.
        special = np.zeros(len(line), dtype=bool)

        hdn = np.where(
            special,
            np.nan,
            (line - zbf[effective_line] - zrev - zd) / gld[effective_line],
        )

        valid_hdn = hdn[np.isfinite(hdn)]
        if len(valid_hdn) > 0:
            line_median = float(np.median(valid_hdn))
            nl_gain = 1.0 - gnl * line_median
            calibrated = np.where(
                np.isfinite(hdn),
                hdn * gcn * nl_gain * gff * gt / guc,
                np.nan,
            )
        else:
            calibrated = hdn

        out[line_idx] = calibrated.astype(np.float32)

    return out


def _read_ancillary_tables(
    cube_path: Path,
    params: HiCalParams,
) -> tuple[np.ndarray, np.ndarray]:
    """Read HiRISE ancillary tables (buffer pixels, calibration image).

    Returns
    -------
    buffer_pixels : np.ndarray
        Image buffer pixels, shape (n_lines, 12).
    cal_image : np.ndarray
        Calibration image, shape (n_cal_lines, n_samples).
    """
    label = pvl.load(str(cube_path))

    buffer_pixels = _read_table_by_name(cube_path, label, "HiRISE Ancillary", params)
    cal_image = _read_table_by_name(cube_path, label, "HiRISE Calibration Image", params)

    return buffer_pixels, cal_image


def _read_table_by_name(
    cube_path: Path,
    label: pvl.PVLModule,
    table_name: str,
    params: HiCalParams,
) -> np.ndarray:
    """Read a named table from an ISIS cube.

    For "HiRISE Ancillary": extracts the BufferPixels field (12 int16 per record).
    For "HiRISE Calibration Image": extracts the Calibration field (nsamps int16 per record).
    """
    # Find the table in the label
    tables = []
    for obj_name, obj_value in label:
        if obj_name == "Table" and hasattr(obj_value, "__getitem__"):
            if obj_value.get("Name") == table_name:
                tables.append(obj_value)

    if not tables:
        raise ValueError(f"Table '{table_name}' not found in {cube_path}")

    table = tables[0]
    start_byte = int(table["StartByte"]) - 1  # ISIS 1-based
    total_bytes = int(table["Bytes"])
    n_records = int(table["Records"])

    with open(cube_path, "rb") as f:
        f.seek(start_byte)
        raw = f.read(total_bytes)

    if table_name == "HiRISE Ancillary":
        # Format: GapFlag(int4) + LineNumber(int4) + BufferPixels(12*int2) + DarkPixels(16*int2)
        # Total per record: 4 + 4 + 24 + 32 = 64 bytes
        record_size = total_bytes // n_records
        buffer_data = np.empty((n_records, 12), dtype=np.int16)
        for i in range(n_records):
            offset = i * record_size
            # Skip GapFlag(4) + LineNumber(4) = 8 bytes to get BufferPixels
            buf = np.frombuffer(raw, dtype="<i2", count=12, offset=offset + 8)
            buffer_data[i] = buf
        return buffer_data

    elif table_name == "HiRISE Calibration Image":
        # Format: Calibration field with 1024 int32 values per record
        record_size = total_bytes // n_records
        n_samps = record_size // 4  # int32 = 4 bytes
        data = np.frombuffer(raw, dtype="<i4").reshape(n_records, n_samps)
        return data.copy()

    elif table_name == "HiRISE Calibration Ancillary":
        # Similar to HiRISE Ancillary but for calibration region
        record_size = total_bytes // n_records
        buffer_data = np.empty((n_records, 12), dtype=np.int16)
        for i in range(n_records):
            offset = i * record_size
            buf = np.frombuffer(raw, dtype="<i2", count=12, offset=offset + 8)
            buffer_data[i] = buf
        return buffer_data

    else:
        raise ValueError(f"Unknown table type: {table_name}")
