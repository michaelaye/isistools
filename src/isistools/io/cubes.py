"""ISIS cube loading and orientation normalization.

Loads ISIS .cub files as xarray DataArrays with proper CRS metadata
and corrected display orientation (north-up for map-projected cubes).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pvl
import rioxarray  # noqa: F401 — registers the .rio accessor
import xarray as xr

# ISIS special pixel constants (IEEE 754 representations).
# Any float32 value <= ISIS_NULL is a special pixel.
ISIS_NULL = np.float32(-3.4028226550889045e38)
ISIS_LRS = np.float32(-3.4028228579130005e38)
ISIS_LIS = np.float32(-3.4028230607370965e38)
ISIS_HRS = np.float32(-3.4028232635611926e38)
ISIS_HIS = np.float32(-3.4028234663852886e38)


def read_label(cube_path: str | Path, *, fast: bool = True) -> pvl.PVLModule:
    """Read the PVL label from an ISIS cube file.

    Parameters
    ----------
    cube_path : path-like
        Path to the .cub file.
    fast : bool
        If True (default), read only the first ``LABEL_MAX_BYTES`` of
        the file and parse that, which is dramatically faster for
        large cubes. ``pvl.load(path)`` without this optimization can
        take seconds on multi-gigabyte files even though the label is
        only a few tens of KB — it apparently does not stop reading
        at the ``End`` marker. For a 2 GB MOLA DEM the difference is
        ~1250 ms vs ~15 ms.

    Returns
    -------
    pvl.PVLModule
        Parsed PVL label.
    """
    if not fast:
        return pvl.load(str(cube_path))

    # Read a generous chunk - ISIS labels are typically < 64 KB but can
    # reach ~1 MB for cubes with many attached Tables or a large History.
    with open(cube_path, "rb") as f:
        head = f.read(LABEL_MAX_BYTES)

    try:
        return pvl.loads(head.decode("latin-1", errors="replace"))
    except Exception:
        # Fall back to the slow path if the chunk truncated the label
        return pvl.load(str(cube_path))


# Upper bound on how many bytes ``read_label(..., fast=True)`` reads.
# Chosen large enough for any ISIS cube we've seen (CTX labels are ~10 KB,
# MOLA DEM ~25 KB, lev1 cubes with Tables ~65 KB). Bumped to 1 MB to cover
# cubes with unusually large History blocks.
LABEL_MAX_BYTES = 1 * 1024 * 1024


def get_projection_info(label: pvl.PVLModule) -> dict | None:
    """Extract map projection information from a cube label.

    Returns None for level-1 (unprojected) cubes.
    """
    try:
        mapping = label["IsisCube"]["Mapping"]
        return dict(mapping)
    except (KeyError, TypeError):
        return None


def get_cube_level(label: pvl.PVLModule) -> int:
    """Determine processing level of an ISIS cube.

    Returns
    -------
    int
        1 for unprojected (raw geometry), 2 for map-projected.
    """
    return 2 if get_projection_info(label) is not None else 1


def load_cube(cube_path: str | Path) -> xr.DataArray:
    """Load an ISIS cube as an xarray DataArray.

    For level-2 (map-projected) cubes, the CRS is set from the
    Mapping group and the image is oriented north-up regardless
    of the original detector readout direction.

    For level-1 cubes, the image is returned in sample/line space
    with label metadata attached as attributes.

    Parameters
    ----------
    cube_path : path-like
        Path to the .cub file.

    Returns
    -------
    xr.DataArray
        Image data with metadata in .attrs.
    """
    cube_path = Path(cube_path)
    if not cube_path.exists():
        raise FileNotFoundError(f"Cube not found: {cube_path}")

    label = read_label(cube_path)

    # rioxarray/rasterio can read ISIS cubes via GDAL's ISIS3 driver
    da = xr.open_dataarray(cube_path, engine="rasterio")

    # Squeeze single-band cubes for simpler downstream handling
    if "band" in da.dims and da.sizes["band"] == 1:
        da = da.squeeze("band", drop=True)

    # Attach useful metadata
    da.attrs["cube_path"] = str(cube_path)
    da.attrs["cube_level"] = get_cube_level(label)
    da.attrs["label"] = label

    proj_info = get_projection_info(label)
    if proj_info is not None:
        da.attrs["projection_info"] = proj_info

    # Normalize orientation for level-2 cubes.
    # ISIS cubes can have flipped line direction depending on the
    # instrument and projection settings. We ensure consistent
    # north-up display by checking the y-coordinate ordering.
    if da.attrs["cube_level"] == 2 and "y" in da.dims:
        y_vals = da.coords["y"].values
        if len(y_vals) > 1 and y_vals[0] < y_vals[-1]:
            # y increases downward — flip to north-up
            da = da.isel(y=slice(None, None, -1))

    return da


def build_serial_lookup(
    cube_paths: list[Path],
) -> dict[str, Path]:
    """Build a mapping from control network serial numbers to cube paths.

    Reads the SpacecraftClockCount from each cube label and matches it
    against the clock count portion of ISIS serial numbers
    (e.g., ``MRO/CTX/0910464726:234``).
    """
    clock_to_path: dict[str, Path] = {}
    for cp in cube_paths:
        try:
            label = pvl.load(str(cp))
            inst = label["IsisCube"]["Instrument"]
            clock = str(
                inst.get("SpacecraftClockCount", inst.get("SpacecraftClockStartCount", ""))
            )
            if clock:
                clock_to_path[clock] = cp
        except Exception:
            continue

    return clock_to_path


def match_serials_to_cubes(
    serial_numbers: list[str],
    cube_paths: list[Path],
) -> dict[str, Path]:
    """Match serial numbers from a control network to cube file paths.

    Returns a dict mapping serial number -> cube path.
    """
    clock_to_path = build_serial_lookup(cube_paths)
    result: dict[str, Path] = {}
    for sn in serial_numbers:
        # Serial numbers are like MRO/CTX/0910464726:234
        # The clock count is the last segment
        clock = sn.rsplit("/", 1)[-1]
        if clock in clock_to_path:
            result[sn] = clock_to_path[clock]
    return result


def get_serial_number(label: pvl.PVLModule) -> str:
    """Construct a serial number proxy from label metadata.

    ISIS serial numbers are typically generated by `serialnumber`
    or stored in the cube's history. This extracts enough info
    to match against control network serial numbers.

    For a proper serial number, use ISIS's `sn` command or the
    SerialNumberList functionality.

    Returns
    -------
    str
        Instrument-derived identifier string.
    """
    inst = label["IsisCube"].get("Instrument", {})
    archive = label["IsisCube"].get("Archive", {})

    # Try common patterns
    spacecraft = inst.get("SpacecraftName", inst.get("SpacecraftId", "Unknown"))
    instrument_id = inst.get("InstrumentId", "Unknown")
    start_time = inst.get("StartTime", inst.get("SpacecraftClockStartCount", "Unknown"))

    return f"{spacecraft}/{instrument_id}/{start_time}"


def _mask_special_pixels(data: np.ndarray) -> None:
    """Replace ISIS special pixel values with NaN, in place."""
    data[data <= np.float32(-3.4028226e38)] = np.nan


def read_isis_cube_raw(cube_path: str | Path) -> tuple[np.ndarray, pvl.PVLModule]:
    """Read raw pixel data and metadata from an ISIS cube.

    Unlike :func:`load_cube` which returns an xarray DataArray via GDAL,
    this reads the binary pixel data directly using numpy. This is used
    by the csm2map processing pipeline where direct array access and
    explicit special-pixel handling are needed.

    Parameters
    ----------
    cube_path : path-like
        Path to an ISIS .cub file.

    Returns
    -------
    data : ndarray, shape (n_lines, n_samples) or (n_bands, n_lines, n_samples)
        Pixel values as float32. ISIS special pixels are converted to NaN.
    label : pvl.PVLModule
        Parsed PVL label.
    """
    cube_path = Path(cube_path)
    label = read_label(cube_path)

    core = label["IsisCube"]["Core"]
    dims = core["Dimensions"]
    n_samples = int(dims["Samples"])
    n_lines = int(dims["Lines"])
    n_bands = int(dims["Bands"])

    fmt = str(core.get("Format", "BandSequential"))

    # For Tile format, use rasterio/GDAL which handles ISIS3 tiling correctly.
    # For BandSequential, reading raw via np.fromfile is faster and avoids the
    # GDAL dependency for that path.
    if fmt.lower() == "tile":
        import rasterio

        with rasterio.open(str(cube_path)) as src:
            data = src.read().astype(np.float32)
            # rasterio returns (bands, rows, cols) which matches our layout
    else:
        pixels = core["Pixels"]
        ptype = str(pixels["Type"]).lower()
        byte_order = str(pixels.get("ByteOrder", "Lsb")).lower()

        dtype_map = {
            "real": np.float32,
            "signedword": np.int16,
            "unsignedbyte": np.uint8,
            "double": np.float64,
        }
        if ptype not in dtype_map:
            msg = f"Unsupported ISIS pixel type: {ptype}"
            raise ValueError(msg)

        np_dtype = dtype_map[ptype]
        endian = "<" if byte_order == "lsb" else ">"
        dt = np.dtype(np_dtype).newbyteorder(endian)

        # ISIS StartByte is 1-based
        start_byte = int(core["StartByte"]) - 1

        data = np.fromfile(
            str(cube_path),
            dtype=dt,
            count=n_bands * n_lines * n_samples,
            offset=start_byte,
        )
        data = data.reshape(n_bands, n_lines, n_samples).astype(np.float32)

    _mask_special_pixels(data)

    if n_bands == 1:
        data = data[0]

    return data, label
