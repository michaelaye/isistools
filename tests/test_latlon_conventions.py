"""Regression tests for the ISIS latitude/longitude convention bug.

Before 0.8.1 the csm2map pipeline read ``LatitudeType``,
``LongitudeDirection`` and ``LongitudeDomain`` from ISIS MAP files but
never acted on them. A MAP file with ``LatitudeType = Planetographic``
was silently treated as planetocentric, producing a ~0.3° latitude
shift on Mars at mid-latitudes — a ~20 km ground-location error with
no warning. The bug was documented as a "latent risk" in the 0.7.0
design doc §6 and then fixed in 0.8.1 by reading the convention
keywords and converting to csm2map's internal planetocentric /
positive-east / 360° convention.

These tests lock in that fix at two levels:

1. **Unit-level**: exercise the pure-geometry helpers in
   ``isistools.geo.projections`` against known Mars values and against
   the sphere identity.

2. **Integration-level**: build two synthetic MAP files that describe
   the same physical ground patch with different conventions
   (Planetocentric vs Planetographic, PositiveEast vs PositiveWest)
   and assert that :func:`grid_from_map_file` produces the same
   csm2map-internal grid from both.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from isistools.geo.projections import (
    normalize_latitude_from_mapping,
    normalize_longitude,
    normalize_longitude_from_mapping,
    planetocentric_to_planetographic,
    planetographic_to_planetocentric,
)
from isistools.processing.grid import grid_from_map_file

# Mars axes (to match F05 and what the Mars MAP files carry)
MARS_EQ = 3396190.0
MARS_POLAR = 3376200.0
# Ratio from the classic conversion: (b/a)^2 ≈ 0.98826
MARS_RATIO_SQ = (MARS_POLAR / MARS_EQ) ** 2


# ------------------------------------------------------------------
# Pure helpers: planetographic ↔ planetocentric


class TestPlanetographicConversion:
    def test_sphere_identity(self):
        """On a spherical body, both conventions are identical."""
        for lat in (-89.0, -45.0, 0.0, 30.0, 90.0):
            assert planetographic_to_planetocentric(lat, 1737.4, 1737.4) == lat
            assert planetocentric_to_planetographic(lat, 1737.4, 1737.4) == lat

    def test_equator_zero_identity(self):
        """At the equator both conventions read 0° exactly."""
        assert planetographic_to_planetocentric(0.0, MARS_EQ, MARS_POLAR) == 0.0
        assert planetocentric_to_planetographic(0.0, MARS_EQ, MARS_POLAR) == 0.0

    def test_poles_unchanged(self):
        """At ±90° both conventions agree (both describe the spin axis)."""
        for lat in (90.0, -90.0):
            assert planetographic_to_planetocentric(lat, MARS_EQ, MARS_POLAR) == lat
            assert planetocentric_to_planetographic(lat, MARS_EQ, MARS_POLAR) == lat

    def test_mars_45_degrees(self):
        """At Mars 45° planetographic, planetocentric is ~44.66° — known value."""
        pc = planetographic_to_planetocentric(45.0, MARS_EQ, MARS_POLAR)
        # Expected: atan((b/a)^2 * tan(45°)) = atan(0.98826) ≈ 44.662°
        expected = math.degrees(math.atan(MARS_RATIO_SQ))
        assert pc == pytest.approx(expected, abs=1e-6)
        assert pc == pytest.approx(44.662, abs=1e-3)

    def test_mars_conversion_is_small_but_nonzero(self):
        """Confirm the magnitude of the Mars-specific shift (the bug's
        numerical impact). The worst-case difference between the two
        conventions on Mars is ~0.33° near ±45°, which at 59 km/° gives
        a ~20 km ground-location error."""
        for pg_lat in (15.0, 30.0, 45.0, 60.0, 75.0):
            pc_lat = planetographic_to_planetocentric(pg_lat, MARS_EQ, MARS_POLAR)
            diff = abs(pg_lat - pc_lat)
            assert 0 < diff < 0.34, f"at {pg_lat}°, diff={diff}°"

    def test_round_trip_scalar(self):
        """pc → pg → pc must be a no-op to floating-point precision."""
        for lat in (-87.3, -60.0, -22.5, 0.0, 17.8, 45.1, 72.4):
            pg = planetocentric_to_planetographic(lat, MARS_EQ, MARS_POLAR)
            pc = planetographic_to_planetocentric(pg, MARS_EQ, MARS_POLAR)
            assert pc == pytest.approx(lat, abs=1e-10)

    def test_round_trip_vectorized(self):
        """Same round-trip test, but vectorized through numpy arrays —
        the ndarray branch of the helper must agree with the scalar
        branch."""
        lats = np.array([-87.3, -60.0, -22.5, 0.0, 17.8, 45.1, 72.4])
        pg = planetocentric_to_planetographic(lats, MARS_EQ, MARS_POLAR)
        pc = planetographic_to_planetocentric(pg, MARS_EQ, MARS_POLAR)
        np.testing.assert_allclose(pc, lats, atol=1e-10)

    def test_vectorized_matches_scalar(self):
        lats = np.array([0.0, 30.0, 45.0, 60.0])
        vec = planetographic_to_planetocentric(lats, MARS_EQ, MARS_POLAR)
        for lat, v in zip(lats.tolist(), vec.tolist()):
            s = planetographic_to_planetocentric(lat, MARS_EQ, MARS_POLAR)
            assert v == pytest.approx(s, abs=1e-12)


# ------------------------------------------------------------------
# Longitude normalization


class TestLongitudeNormalization:
    def test_positive_east_360_is_noop(self):
        assert normalize_longitude(120.0, direction="PositiveEast", domain=360) == 120.0
        assert normalize_longitude(0.0, direction="PositiveEast", domain=360) == 0.0

    def test_positive_east_180_wraps(self):
        """A value in [-180, 180] should wrap into [0, 360) when we
        normalize to the 360 domain."""
        assert normalize_longitude(-10.0, direction="PositiveEast", domain=180) == pytest.approx(
            350.0
        )
        assert normalize_longitude(170.0, direction="PositiveEast", domain=180) == pytest.approx(
            170.0
        )

    def test_positive_west_flips_sign(self):
        """PositiveWest → PositiveEast: lon_pe = -lon_pw (mod 360)."""
        assert normalize_longitude(287.05, direction="PositiveWest", domain=360) == pytest.approx(
            72.95, abs=1e-6
        )
        assert normalize_longitude(72.95, direction="PositiveWest", domain=360) == pytest.approx(
            287.05, abs=1e-6
        )

    def test_positive_west_round_trip(self):
        for lon in (0.0, 30.0, 90.0, 180.0, 270.0, 359.99):
            pe = normalize_longitude(lon, direction="PositiveWest", domain=360)
            back = normalize_longitude(pe, direction="PositiveWest", domain=360)
            assert back == pytest.approx(lon, abs=1e-9)

    def test_unknown_direction_raises(self):
        with pytest.raises(ValueError, match="Unrecognized longitude direction"):
            normalize_longitude(0.0, direction="Sideways", domain=360)


# ------------------------------------------------------------------
# Mapping-aware wrappers


class TestMappingHelpers:
    def test_normalize_latitude_default_is_planetocentric(self):
        """Mapping with no LatitudeType is assumed to be Planetocentric
        (ISIS default)."""
        mapping = {}
        assert normalize_latitude_from_mapping(45.0, mapping, MARS_EQ, MARS_POLAR) == 45.0

    def test_normalize_latitude_planetographic_converts(self):
        mapping = {"LatitudeType": "Planetographic"}
        pc = normalize_latitude_from_mapping(45.0, mapping, MARS_EQ, MARS_POLAR)
        assert pc == pytest.approx(44.662, abs=1e-3)

    def test_normalize_latitude_explicit_planetocentric_noop(self):
        mapping = {"LatitudeType": "Planetocentric"}
        assert normalize_latitude_from_mapping(30.0, mapping, MARS_EQ, MARS_POLAR) == 30.0

    def test_normalize_latitude_rejects_junk(self):
        mapping = {"LatitudeType": "Flatland"}
        with pytest.raises(ValueError, match="Unrecognized LatitudeType"):
            normalize_latitude_from_mapping(0.0, mapping, MARS_EQ, MARS_POLAR)

    def test_normalize_longitude_default_is_positive_east(self):
        mapping = {}
        assert normalize_longitude_from_mapping(120.0, mapping) == 120.0

    def test_normalize_longitude_positive_west(self):
        mapping = {"LongitudeDirection": "PositiveWest", "LongitudeDomain": 360}
        assert normalize_longitude_from_mapping(287.05, mapping) == pytest.approx(72.95, abs=1e-6)


# ------------------------------------------------------------------
# Integration: grid_from_map_file honors the conventions


def _write_map_file(
    path: Path,
    *,
    lat_type: str,
    lon_direction: str,
    lon_domain: int,
    min_lat: float,
    max_lat: float,
    min_lon: float,
    max_lon: float,
    resolution: float = 6.0,
) -> None:
    """Write a minimal ISIS MAP PVL file with the given conventions."""
    content = f"""Group = Mapping
  ProjectionName     = Equirectangular
  CenterLatitude     = 0.0
  CenterLongitude    = 0.0
  TargetName         = Mars
  EquatorialRadius   = {MARS_EQ:.1f} <meters>
  PolarRadius        = {MARS_POLAR:.1f} <meters>
  LatitudeType       = {lat_type}
  LongitudeDirection = {lon_direction}
  LongitudeDomain    = {lon_domain}
  MinimumLatitude    = {min_lat}
  MaximumLatitude    = {max_lat}
  MinimumLongitude   = {min_lon}
  MaximumLongitude   = {max_lon}
  PixelResolution    = {resolution} <meters/pixel>
End_Group
"""
    path.write_text(content)


class TestGridFromMapFileConventions:
    """Integration test: the same physical ground patch described in two
    different conventions must produce the same csm2map-internal grid
    (planetocentric, positive-east, 360°)."""

    def test_planetocentric_vs_planetographic_produce_same_grid(self, tmp_path):
        """Build one MAP file in each convention describing the SAME
        physical patch. After conversion, csm2map should see identical
        planetocentric lat ranges from both.

        Before the 0.8.1 fix this test would fail: the planetographic
        file's lat values would be used verbatim as planetocentric,
        producing a grid shifted by ~0.3° at mid-latitudes."""
        # Planetocentric reference values for a mid-latitude Mars patch.
        pc_min = 30.0
        pc_max = 45.0
        # Convert to planetographic (what an ISIS MAP file with
        # LatitudeType=Planetographic would carry for the same patch).
        pg_min = planetocentric_to_planetographic(pc_min, MARS_EQ, MARS_POLAR)
        pg_max = planetocentric_to_planetographic(pc_max, MARS_EQ, MARS_POLAR)

        map_pc = tmp_path / "planetocentric.map"
        map_pg = tmp_path / "planetographic.map"

        _write_map_file(
            map_pc,
            lat_type="Planetocentric",
            lon_direction="PositiveEast",
            lon_domain=360,
            min_lat=pc_min,
            max_lat=pc_max,
            min_lon=100.0,
            max_lon=110.0,
        )
        _write_map_file(
            map_pg,
            lat_type="Planetographic",
            lon_direction="PositiveEast",
            lon_domain=360,
            min_lat=pg_min,
            max_lat=pg_max,
            min_lon=100.0,
            max_lon=110.0,
        )

        grid_pc = grid_from_map_file(map_pc)

        # The planetographic MAP file emits a conversion warning we
        # want to see — it confirms the fix's user-facing signal.
        with pytest.warns(UserWarning, match="LatitudeType='Planetographic'"):
            grid_pg = grid_from_map_file(map_pg)

        # After conversion, the internal lat_min / lat_max must match.
        assert grid_pg.lat_min == pytest.approx(grid_pc.lat_min, abs=1e-9)
        assert grid_pg.lat_max == pytest.approx(grid_pc.lat_max, abs=1e-9)
        # And the grid dimensions must match to the pixel — if we were
        # off by the ~0.33° Mars planetographic/centric shift, the row
        # count would differ by ~3% at these latitudes.
        assert grid_pg.height == grid_pc.height
        assert grid_pg.width == grid_pc.width

    def test_positive_west_vs_positive_east_produce_same_grid(self, tmp_path):
        """Longitude convention flip: PositiveWest → PositiveEast must
        produce the same physical patch."""
        # Positive-east reference values
        pe_min = 100.0
        pe_max = 110.0
        # Positive-west values for the same patch: lon_pw = 360 - lon_pe,
        # with min/max swapping.
        pw_min = 360.0 - pe_max
        pw_max = 360.0 - pe_min

        map_pe = tmp_path / "positive_east.map"
        map_pw = tmp_path / "positive_west.map"

        _write_map_file(
            map_pe,
            lat_type="Planetocentric",
            lon_direction="PositiveEast",
            lon_domain=360,
            min_lat=30.0,
            max_lat=45.0,
            min_lon=pe_min,
            max_lon=pe_max,
        )
        _write_map_file(
            map_pw,
            lat_type="Planetocentric",
            lon_direction="PositiveWest",
            lon_domain=360,
            min_lat=30.0,
            max_lat=45.0,
            min_lon=pw_min,
            max_lon=pw_max,
        )

        grid_pe = grid_from_map_file(map_pe)
        with pytest.warns(UserWarning, match="LongitudeDirection='PositiveWest'"):
            grid_pw = grid_from_map_file(map_pw)

        assert grid_pw.lon_min == pytest.approx(grid_pe.lon_min, abs=1e-9)
        assert grid_pw.lon_max == pytest.approx(grid_pe.lon_max, abs=1e-9)
        assert grid_pw.width == grid_pe.width
        assert grid_pw.height == grid_pe.height

    def test_planetocentric_positive_east_silent(self, tmp_path):
        """The common case (what csm2map writes itself) must NOT emit
        any conversion warning — the convention matches the default."""
        import warnings

        map_path = tmp_path / "default.map"
        _write_map_file(
            map_path,
            lat_type="Planetocentric",
            lon_direction="PositiveEast",
            lon_domain=360,
            min_lat=30.0,
            max_lat=45.0,
            min_lon=100.0,
            max_lon=110.0,
        )

        with warnings.catch_warnings():
            warnings.simplefilter("error")  # turn any UserWarning into an exception
            # ... except for unrelated pyproj warnings, which are okay
            warnings.filterwarnings("ignore", category=UserWarning, module="pyproj")
            grid = grid_from_map_file(map_path)

        assert grid.lat_min == pytest.approx(30.0, abs=1e-9)
        assert grid.lat_max == pytest.approx(45.0, abs=1e-9)
