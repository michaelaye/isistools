"""DEM-based surface radius sampling for the csm2map pipeline.

ISIS uses an ellipsoidal target model by default but allows a Digital
Elevation Model (DEM) to be supplied as the shape model via the cube
label's ``Kernels.ShapeModel`` keyword. The default Mars DEM is the MOLA
``base/dems/molaMarsPlanetaryRadius0005.cub`` which stores per-pixel
local body radius (not elevation) as ``int16`` with ``Base=3396000``.

This module bilinearly samples that DEM at arbitrary lat/lon arrays and
returns the local radius in meters. The sampler reads only the window of
the DEM that covers the requested points (so a small CTX image only
loads a few hundred kilobytes from a 2 GB DEM file), making it cheap to
plug into ``compute_transform_coarse`` / ``compute_transform_dense``
where we need a per-point radius for ``ground_to_image_batch``.
"""

import os
from pathlib import Path

import numpy as np
import rasterio
from pyproj import CRS, Transformer
from scipy.ndimage import map_coordinates

from isistools.csm2map.projections import _to_meters, mapping_to_crs
from isistools.io.cubes import read_label

# pyproj raises when transforming between Earth (EPSG:4326) and a Mars CRS;
# we have to opt out of that check explicitly because both sides describe
# planetocentric lat/lon for the same body, just with different axes.
os.environ.setdefault("PROJ_IGNORE_CELESTIAL_BODY", "YES")


class DemRadiusSampler:
    """Bilinearly sample local body radius from an ISIS DEM cube.

    The DEM must be an ISIS cube whose ``Mapping`` group describes a
    standard equirectangular / sinusoidal / simple-cylindrical projection
    and whose pixel values represent body radius in meters (after
    applying ``Base`` and ``Multiplier``).

    Parameters
    ----------
    dem_path : path-like
        Path to the DEM cube (e.g. ``$ISISDATA/base/dems/molaMarsPlanetaryRadius0005.cub``).
    fallback_radius : float
        Radius (m) to use for points where the DEM has nodata or for points
        outside the DEM coverage. **Required, no default.** Callers should
        pass the target body's mean radius (``TargetBody.radius_mean_m``);
        hardcoding a Mars-specific default would silently miss-project
        non-Mars cubes.

    Notes
    -----
    The sampler keeps the DEM file open between calls (rasterio handle)
    and caches a single windowed read of the most-recently-requested
    bounding box. Repeated calls within the same lat/lon region are O(1)
    file IO; switching to a different region triggers a new window read.
    """

    def __init__(self, dem_path: str | Path, fallback_radius: float) -> None:
        self.dem_path = Path(dem_path)
        if not self.dem_path.exists():
            msg = f"DEM cube not found: {self.dem_path}"
            raise FileNotFoundError(msg)

        # Parse the label to get pixel scaling and projection metadata
        label = read_label(self.dem_path)
        core = label["IsisCube"]["Core"]
        mapping = label["IsisCube"]["Mapping"]
        pixels = core["Pixels"]

        self.base = float(pixels.get("Base", 0.0))
        self.mult = float(pixels.get("Multiplier", 1.0))
        self.nodata = -32768  # ISIS SignedWord NULL — also rasterio nodata

        dims = core["Dimensions"]
        self.n_samples = int(dims["Samples"])
        self.n_lines = int(dims["Lines"])

        self.crs = mapping_to_crs(mapping)
        self.x_ul = _to_meters(mapping["UpperLeftCornerX"])
        self.y_ul = _to_meters(mapping["UpperLeftCornerY"])
        self.pixel_resolution = _to_meters(mapping["PixelResolution"])

        self.fallback_radius = float(fallback_radius)

        # rasterio handle (lazy open)
        self._src = None
        # Cached window
        self._win_data: np.ndarray | None = None  # float32 radii
        self._win_col_ofs = 0
        self._win_row_ofs = 0
        self._win_width = 0
        self._win_height = 0

        # Reusable transformer
        self._tr_to_dem = Transformer.from_crs(
            CRS.from_epsg(4326),
            self.crs,
            always_xy=True,
        )

    def __del__(self) -> None:
        if self._src is not None:
            try:
                self._src.close()
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Coordinate helpers

    def _latlon_to_pixel(
        self,
        lats_deg: np.ndarray,
        lons_deg: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Convert lat/lon (degrees) to fractional DEM (col, row)."""
        x_map, y_map = self._tr_to_dem.transform(lons_deg, lats_deg)
        col = (x_map - self.x_ul) / self.pixel_resolution - 0.5
        row = (self.y_ul - y_map) / self.pixel_resolution - 0.5
        return col, row

    # ------------------------------------------------------------------
    # Window IO

    def _ensure_window(
        self,
        col_min: float,
        col_max: float,
        row_min: float,
        row_max: float,
    ) -> None:
        """Read the DEM window covering the requested bbox + buffer."""
        # Buffer for bilinear kernel + safety
        buf = 4
        col0 = max(0, int(np.floor(col_min)) - buf)
        col1 = min(self.n_samples, int(np.ceil(col_max)) + buf + 1)
        row0 = max(0, int(np.floor(row_min)) - buf)
        row1 = min(self.n_lines, int(np.ceil(row_max)) + buf + 1)

        if col1 <= col0 or row1 <= row0:
            # Whole bbox is outside the DEM
            self._win_data = None
            return

        # If the cached window already covers this region, reuse it
        if (
            self._win_data is not None
            and col0 >= self._win_col_ofs
            and col1 <= self._win_col_ofs + self._win_width
            and row0 >= self._win_row_ofs
            and row1 <= self._win_row_ofs + self._win_height
        ):
            return

        if self._src is None:
            self._src = rasterio.open(str(self.dem_path))

        win = rasterio.windows.Window(col0, row0, col1 - col0, row1 - row0)
        raw = self._src.read(1, window=win)  # int16
        # Convert to radii in meters; mark nodata as NaN
        radii = raw.astype(np.float32) * self.mult + self.base
        radii[raw == self.nodata] = np.nan

        self._win_data = radii
        self._win_col_ofs = col0
        self._win_row_ofs = row0
        self._win_width = col1 - col0
        self._win_height = row1 - row0

    # ------------------------------------------------------------------
    # Public API

    def sample_radii(
        self,
        lats_rad: np.ndarray,
        lons_rad: np.ndarray,
    ) -> np.ndarray:
        """Bilinearly sample local body radius (m) at the given points.

        Parameters
        ----------
        lats_rad, lons_rad : ndarray
            Planetocentric latitude and longitude in radians. Any shape;
            output shape matches the inputs.

        Returns
        -------
        radii : ndarray of float
            Local body radius in meters at each input point. Points that
            fall outside the DEM or hit nodata get ``fallback_radius``.
        """
        in_shape = lats_rad.shape
        lats_deg = np.rad2deg(lats_rad).ravel()
        lons_deg = np.rad2deg(lons_rad).ravel()

        col, row = self._latlon_to_pixel(lats_deg, lons_deg)

        col_min = float(np.nanmin(col)) if col.size else 0.0
        col_max = float(np.nanmax(col)) if col.size else 0.0
        row_min = float(np.nanmin(row)) if row.size else 0.0
        row_max = float(np.nanmax(row)) if row.size else 0.0

        self._ensure_window(col_min, col_max, row_min, row_max)

        if self._win_data is None:
            return np.full(in_shape, self.fallback_radius, dtype=np.float64)

        # Local indices into the windowed array
        local_col = col - self._win_col_ofs
        local_row = row - self._win_row_ofs

        coords = np.array([local_row, local_col])
        radii = map_coordinates(
            self._win_data.astype(np.float64),
            coords,
            order=1,  # bilinear; cubic risks ringing on flat DEM
            mode="constant",
            cval=np.nan,
        )

        # Replace any NaN (from nodata or out-of-window) with fallback
        bad = ~np.isfinite(radii)
        if bad.any():
            radii[bad] = self.fallback_radius

        return radii.reshape(in_shape)


def resolve_shape_model(cube_path: str | Path) -> Path | None:
    """Resolve a cube's ``Kernels.ShapeModel`` keyword to a real DEM path.

    Reads the input cube's PVL label, extracts the ``ShapeModel`` keyword
    if present, and expands ``$base``/``$ISISDATA`` to the actual on-disk
    path using the ``ISISDATA`` environment variable.

    Returns
    -------
    Path or None
        The resolved DEM cube path if the cube specifies a real DEM
        (e.g. ``$base/dems/molaMarsPlanetaryRadius0005.cub``), otherwise
        ``None`` (e.g. for cubes with ``ShapeModel = Null`` or ``Ellipsoid``).
    """
    label = read_label(cube_path)
    kernels = label["IsisCube"].get("Kernels", {})
    raw = kernels.get("ShapeModel")
    if raw is None:
        return None

    text = str(raw).strip()
    if not text or text.lower() in ("null", "none", "ellipsoid", "system"):
        return None

    # Expand $ISISDATA / $base
    isisdata = os.environ.get("ISISDATA")
    if isisdata:
        text = text.replace("$ISISDATA", isisdata).replace("$base", f"{isisdata}/base")

    path = Path(text)
    if not path.exists():
        return None
    return path
