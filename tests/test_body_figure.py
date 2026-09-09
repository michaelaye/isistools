"""Regression tests for the ocentric-latitudes-on-the-ographic-ellipsoid bug.

Through 0.11.0, every CRS csm2map built emitted ``+a=EquatorialRadius
+b=PolarRadius`` — the *ographic* ellipsoid — while the latitudes being
projected were planetocentric (from the camera model, or converted to
planetocentric by ``grid_from_map_file``). PROJ then applied the ellipsoidal
meridional scaling to latitudes that were never ographic.

On Mars that is a fixed ``1 - e² = -1.174 %`` along-track error in any
projection that honours flattening. Measured against ISIS ``cam2map`` over five
CTX cubes spanning 7 168 to 52 224 lines: **-1.1709 %, sd 0.0107** in
Sinusoidal, reaching 3.7 km on a full-length strip. Equirectangular hid it
entirely, because PROJ's ``eqc`` ignores flattening — the grid was right and
only the datum label was wrong.

It is also a CRS error, not just a scale one. Planetocentric Mars is
``IAU_2015:49900`` (*Mars (2015) - Sphere / Ocentric*, ``f = 0``); ``+a``/``+b``
describes ``49901``, the ographic ellipsoid.

These tests lock in:

1. ``body_figure_params`` picking the figure from the latitude convention.
2. ``mapping_to_crs`` putting a planetocentric MAP file on the sphere even when
   the Mapping group carries a genuine ``PolarRadius``.
3. The numeric consequence in Sinusoidal — the case the -1.17 % was measured in.
4. Equirectangular staying byte-identical, so the fix is not a silent re-grid of
   existing equirectangular products.
"""

from __future__ import annotations

import math

import pytest
from pyproj import CRS, Transformer

from isistools.csm2map.projections import (
    body_figure_params,
    is_planetocentric,
    mapping_to_crs,
)

# IAU 2015 Mars
A = 3396190.0
B = 3376200.0
E2 = (A**2 - B**2) / A**2  # 0.0117374


def _mapping(proj_name: str, lat_type: str | None = "Planetocentric") -> dict:
    m = {
        "ProjectionName": proj_name,
        "EquatorialRadius": A,
        "PolarRadius": B,
        "CenterLongitude": 0.0,
        "CenterLatitude": 0.0,
    }
    if lat_type is not None:
        m["LatitudeType"] = lat_type
    return m


class TestBodyFigureParams:
    @pytest.mark.parametrize(
        "lat_type", ["Planetocentric", "planetocentric", "centric", "  Centric  "]
    )
    def test_centric_gets_the_sphere(self, lat_type):
        assert body_figure_params(A, B, lat_type) == [f"+R={A}"]

    def test_absent_latitude_type_defaults_to_centric(self):
        """Matches grid_from_map_file's default, which is also Planetocentric."""
        assert body_figure_params(A, B, None) == [f"+R={A}"]
        assert body_figure_params(A, B) == [f"+R={A}"]

    @pytest.mark.parametrize("lat_type", ["Planetographic", "planetographic", "ographic"])
    def test_ographic_keeps_the_ellipsoid(self, lat_type):
        assert body_figure_params(A, B, lat_type) == [f"+a={A}", f"+b={B}"]

    def test_is_planetocentric(self):
        assert is_planetocentric("Planetocentric")
        assert is_planetocentric(None)
        assert not is_planetocentric("Planetographic")


class TestIsisIsNotUniformlySpherical:
    """ISIS computes some projections on a sphere and others on the ellipsoid.

    Measured by projecting ``B04_011267_0983_XN_81S063W`` (81-82 deg S) with a
    **planetocentric** PolarStereographic MAP file and probing with ``mappt`` at
    planetocentric -81.8 deg, lon 296. ISIS returned
    ``X = -381411.695, Y = -297991.475`` (radius 484018.388 m), which matches
    ``ellipsoid x ographic latitude`` to the millimetre. The alternatives:
    sphere x ocentric is +0.5920 %, sphere x ographic -0.5768 %, ellipsoid x
    ocentric +1.1753 %.

    Getting this wrong costs ~2.9 km at 81.8 deg S, and PolarStereographic is
    what ``ctxpipe._auto_projection`` selects above |lat| 70.
    """

    ISIS_RADIUS_M = 484018.388
    LAT_OC, LAT_OG, LON = -81.8, -81.894960151079, 296.0

    def test_polar_stereographic_keeps_the_ellipsoid(self):
        assert body_figure_params(A, B, "Planetocentric", "PolarStereographic") == [
            f"+a={A}",
            f"+b={B}",
        ]

    def test_spherical_projections_still_get_the_sphere(self):
        for proj in ("Equirectangular", "SimpleCylindrical", "Sinusoidal", "Mercator"):
            assert body_figure_params(A, B, "Planetocentric", proj) == [f"+R={A}"]

    def test_unverified_projection_keeps_historical_ellipsoid(self):
        """Never widen the spherical set without measuring against ISIS first."""
        for proj in ("LambertConformal", "Orthographic", "Robinson"):
            assert body_figure_params(A, B, "Planetocentric", proj) == [f"+a={A}", f"+b={B}"]

    def test_reproduces_the_isis_polar_measurement(self):
        figure = " ".join(body_figure_params(A, B, "Planetocentric", "PolarStereographic"))
        crs = CRS.from_proj4(
            f"+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=64 {figure} +units=m +no_defs +type=crs"
        )
        geo = CRS.from_proj4(f"+proj=longlat {figure} +no_defs +type=crs")
        x, y = Transformer.from_crs(geo, crs, always_xy=True).transform(self.LON, self.LAT_OG)
        assert math.hypot(x, y) == pytest.approx(self.ISIS_RADIUS_M, abs=0.01)

    def test_sphere_would_miss_the_isis_polar_measurement(self):
        """Guard against reverting to a blanket sphere rule."""
        crs = CRS.from_proj4(
            f"+proj=stere +lat_0=-90 +lat_ts=-90 +lon_0=64 +R={A} +units=m +no_defs +type=crs"
        )
        geo = CRS.from_proj4(f"+proj=longlat +R={A} +no_defs +type=crs")
        x, y = Transformer.from_crs(geo, crs, always_xy=True).transform(self.LON, self.LAT_OC)
        assert abs(math.hypot(x, y) - self.ISIS_RADIUS_M) > 2000.0


class TestMappingToCrs:
    def test_planetocentric_mapping_is_spherical(self):
        """The bug: a real PolarRadius in the Mapping group used to win."""
        crs = mapping_to_crs(_mapping("Equirectangular"))
        assert crs.ellipsoid.semi_major_metre == pytest.approx(A)
        assert crs.ellipsoid.semi_minor_metre == pytest.approx(A)
        assert crs.ellipsoid.inverse_flattening in (0.0, math.inf, None)

    def test_planetocentric_matches_iau_49900_figure(self):
        """The registered CRS for planetocentric Mars is the 49900 sphere."""
        ours = mapping_to_crs(_mapping("Equirectangular")).ellipsoid
        iau = CRS.from_user_input("IAU_2015:49900").ellipsoid
        assert ours.semi_major_metre == pytest.approx(iau.semi_major_metre)
        assert ours.semi_minor_metre == pytest.approx(iau.semi_minor_metre)

    def test_planetographic_mapping_keeps_the_ellipsoid(self):
        crs = mapping_to_crs(_mapping("Equirectangular", "Planetographic"))
        assert crs.ellipsoid.semi_major_metre == pytest.approx(A)
        assert crs.ellipsoid.semi_minor_metre == pytest.approx(B)


class TestGeometricConsequence:
    """The -1.17 %: only visible where PROJ honours flattening."""

    LAT = 0.7733354  # the G09_021641 latitude span used in the measurement

    def _y(self, crs: CRS, lat: float) -> float:
        tf = Transformer.from_crs(CRS.from_epsg(4326), crs, always_xy=True)
        return tf.transform(0.0, lat)[1]

    def test_sinusoidal_uses_the_spherical_scale(self):
        y = self._y(mapping_to_crs(_mapping("Sinusoidal")), self.LAT)
        expected_sphere = A * math.radians(self.LAT)  # 45839.2 m — matches cam2map
        assert y == pytest.approx(expected_sphere, rel=1e-4)

    def test_sinusoidal_no_longer_matches_the_ellipsoidal_scale(self):
        """Guard against a silent revert: the two differ by 1 - e²."""
        y = self._y(mapping_to_crs(_mapping("Sinusoidal")), self.LAT)
        ellipsoidal = A * (1 - E2) * math.radians(self.LAT)  # ~45301 m — the bug
        assert abs(y - ellipsoidal) / ellipsoidal > 0.01
        assert (ellipsoidal - y) / y == pytest.approx(-E2, abs=2e-4)

    def test_equirectangular_pins_to_the_spherical_scale(self):
        """Equirectangular is figure-sensitive too, as of PROJ 9.8.

        Measured in this environment: PROJ **9.7.1** returned ``A*rad(lat)`` for
        both ``+R`` and ``+a/+b`` — ``eqc`` was spherical-only, which is why
        equirectangular was previously described as immune to the figure choice.
        PROJ **9.8.1** returns the true meridian arc for ``+a/+b``, a difference
        of ``1 - e²`` (-1.174 % on Mars).

        So declaring an ellipsoid you did not mean now shifts equirectangular
        geometry *as a function of the PROJ version* — silently, on a library
        upgrade. Pinning the planetocentric path to the sphere makes the output
        version-independent, which is the guarantee this fix actually provides,
        and is asserted here rather than any sphere-vs-ellipsoid equality.
        """
        centric = self._y(mapping_to_crs(_mapping("Equirectangular")), self.LAT)
        assert centric == pytest.approx(A * math.radians(self.LAT), rel=1e-9)


class TestMappingPvlSidecar:
    """The sidecar must describe the same figure of the body as the GeoTIFF.

    ``write_mapping_pvl`` used to emit the *kernel* radii (``body.radius_*_m``)
    regardless of the CRS the data was projected with. After the figure fix, a
    planetocentric Equirectangular product sits on the 49900 sphere, so the
    sidecar was claiming ``PolarRadius = 3376200`` — 49901's ographic
    semi-minor, 20 km from the file's own datum.

    That mattered most exactly where the sidecar is most likely to be consulted:
    when a reader has degraded the GeoTIFF's CRS to ``None`` (WKT1 has no slot
    for the IAU edition) and someone reaches for the plain-text file to
    disambiguate 49900 from 49901.
    """

    @staticmethod
    def _body():
        from isistools.csm2map.camera import TargetBody

        return TargetBody(
            name="MARS",
            naif_id=499,
            radius_equatorial_m=A,
            radius_polar_m=B,  # the kernel ellipsoid the writer used to emit
            radius_mean_m=(2.0 * A + B) / 3.0,
        )

    @staticmethod
    def _grid(crs):
        from isistools.csm2map.grid import grid_from_params

        return grid_from_params(
            crs=crs, resolution=5.0,
            lat_min=-1.0, lat_max=1.0, lon_min=-1.0, lon_max=1.0,
        )

    @staticmethod
    def _read(pvl_path):
        out = {}
        for line in pvl_path.read_text().splitlines():
            if "=" in line:
                k, _, v = line.partition("=")
                out[k.strip()] = v.strip().split("<")[0].strip()
        return out

    def test_spherical_crs_yields_spherical_radii(self, tmp_path):
        from isistools.csm2map.writers import write_mapping_pvl

        crs = CRS.from_user_input("IAU_2015:49910")  # sphere / ocentric / eqc
        p = write_mapping_pvl(tmp_path / "x.tif", self._grid(crs), self._body())
        got = self._read(p)
        assert float(got["EquatorialRadius"]) == pytest.approx(A)
        # The regression: this used to be B (3376200), the ographic semi-minor.
        assert float(got["PolarRadius"]) == pytest.approx(A)

    def test_sidecar_agrees_with_the_crs_it_sits_beside(self, tmp_path):
        """Holds for the ellipsoidal case too — the sidecar follows the CRS."""
        from isistools.csm2map.writers import write_mapping_pvl

        for crs in (CRS.from_user_input("IAU_2015:49910"),
                    CRS.from_proj4(f"+proj=eqc +a={A} +b={B} +no_defs +type=crs")):
            p = write_mapping_pvl(tmp_path / "y.tif", self._grid(crs), self._body())
            got = self._read(p)
            assert float(got["EquatorialRadius"]) == pytest.approx(
                crs.ellipsoid.semi_major_metre)
            assert float(got["PolarRadius"]) == pytest.approx(
                crs.ellipsoid.semi_minor_metre)
