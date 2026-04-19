"""Regression tests for the body-agnostic csm2map refactor.

Before 0.8.0 the csm2map pipeline hardcoded Mars radii in four places
(camera.py::get_target_radii, project.py's default projection string,
dem.py's fallback_radius default, and projections.py::mapping_to_crs).
A JANUS, LRO, or Europa Clipper cube passed to csm2map would silently
project with Mars radii → ~6× scale error with no warning.

This file locks in the regression fixes by exercising:

1. :class:`TargetBody` correctness across multiple bodies — Mars, Moon,
   Europa, and a hypothetical "BODY999" that no lookup table could have
   known about. The body info always flows from the ALE ISD, never from
   a body-name allowlist.

2. Unit conversion (km → m) from the ISD ``radii`` block.

3. Cross-check between ``isd["radii"]`` and
   ``isd["naif_keywords"]["BODY<code>_RADII"]`` — the two should agree
   to meter precision because they come from the same SPICE kernel load.

4. :func:`mapping_to_crs` fails loudly when given a Mapping group with
   no explicit radii — it previously silently defaulted to Mars, which
   is the exact failure mode this refactor eliminates.
"""

from __future__ import annotations

import pytest

from isistools.csm2map.camera import TargetBody
from isistools.csm2map.projections import mapping_to_crs


def _make_isd(naif_id: int, semimajor_km: float, semiminor_km: float) -> dict:
    """Build a minimal ALE-shaped ISD dict for a given body."""
    # semimajor comes first and third, semiminor second in the 3-vector
    # (SPICE BODY_RADII convention: equatorial_x, equatorial_y, polar)
    body_radii_km = [semimajor_km, semimajor_km, semiminor_km]
    return {
        "radii": {
            "semimajor": semimajor_km,
            "semiminor": semiminor_km,
            "unit": "km",
        },
        "naif_keywords": {
            "BODY_CODE": naif_id,
            f"BODY{naif_id}_RADII": body_radii_km,
        },
    }


# ------------------------------------------------------------------
# TargetBody.from_isd: per-body correctness


class TestTargetBodyFromIsd:
    def test_mars(self):
        """Mars: the body we had hardcoded — must produce the exact same numbers."""
        isd = _make_isd(499, semimajor_km=3396.19, semiminor_km=3376.2)
        body = TargetBody.from_isd(isd, target_name="Mars")

        assert body.name == "MARS"
        assert body.naif_id == 499
        assert body.radius_equatorial_m == pytest.approx(3396190.0, abs=1.0)
        assert body.radius_polar_m == pytest.approx(3376200.0, abs=1.0)
        # Mean radius: (2*a + c) / 3 = (2*3396190 + 3376200) / 3 ≈ 3389526.67
        assert body.radius_mean_m == pytest.approx(3389526.67, abs=1.0)

    def test_moon(self):
        """Moon: different body, different radii — must be correctly extracted."""
        isd = _make_isd(301, semimajor_km=1737.4, semiminor_km=1737.4)
        body = TargetBody.from_isd(isd, target_name="Moon")

        assert body.name == "MOON"
        assert body.naif_id == 301
        assert body.radius_equatorial_m == pytest.approx(1737400.0)
        assert body.radius_polar_m == pytest.approx(1737400.0)
        # Spherical body: mean == equatorial == polar
        assert body.radius_mean_m == pytest.approx(1737400.0)

    def test_europa(self):
        """Europa: JUICE/Clipper target — must NOT leak Mars radii."""
        isd = _make_isd(502, semimajor_km=1562.6, semiminor_km=1559.5)
        body = TargetBody.from_isd(isd, target_name="Europa")

        assert body.name == "EUROPA"
        assert body.naif_id == 502
        assert body.radius_equatorial_m == pytest.approx(1562600.0)
        assert body.radius_polar_m == pytest.approx(1559500.0)
        # Critical: mean radius must NOT be anywhere near Mars's 3.39 Mm.
        assert body.radius_mean_m < 1_600_000.0
        assert body.radius_mean_m > 1_500_000.0

    def test_unknown_body_synthesizes_name(self):
        """A hypothetical new target with no target_name provided — the body
        must still work. This is the test that would have caught the old
        ``BODY499_RADII`` hardcode: a body with NAIF ID 999 that no lookup
        table could know about still produces correct radii."""
        isd = _make_isd(999, semimajor_km=100.0, semiminor_km=80.0)
        body = TargetBody.from_isd(isd)  # no target_name

        assert body.name == "BODY999"
        assert body.naif_id == 999
        assert body.radius_equatorial_m == pytest.approx(100_000.0)
        assert body.radius_polar_m == pytest.approx(80_000.0)

    def test_target_name_uppercased(self):
        """Lowercase or mixed-case target names from the cube label get
        normalized to uppercase to match SPICE conventions."""
        isd = _make_isd(499, 3396.19, 3376.2)
        body = TargetBody.from_isd(isd, target_name="mars")
        assert body.name == "MARS"


# ------------------------------------------------------------------
# TargetBody.from_isd: error handling


class TestTargetBodyFromIsdErrors:
    def test_unit_meters_supported(self):
        """Some ISD variants might report radii in meters directly. Don't
        corrupt them by re-multiplying by 1000."""
        isd = {
            "radii": {
                "semimajor": 3396190.0,
                "semiminor": 3376200.0,
                "unit": "m",
            },
            "naif_keywords": {
                "BODY_CODE": 499,
                # BODY499_RADII also in meters when unit=='m' — the
                # cross-check multiplies by the same scale.
                "BODY499_RADII": [3396190.0, 3396190.0, 3376200.0],
            },
        }
        body = TargetBody.from_isd(isd, target_name="Mars")
        assert body.radius_equatorial_m == pytest.approx(3396190.0, abs=1.0)

    def test_unknown_unit_raises(self):
        isd = {
            "radii": {"semimajor": 3396.19, "semiminor": 3376.2, "unit": "furlongs"},
            "naif_keywords": {"BODY_CODE": 499, "BODY499_RADII": [3396.19, 3396.19, 3376.2]},
        }
        with pytest.raises(ValueError, match="Unrecognized radii unit"):
            TargetBody.from_isd(isd)

    def test_missing_body_code_raises(self):
        isd = {
            "radii": {"semimajor": 3396.19, "semiminor": 3376.2, "unit": "km"},
            "naif_keywords": {},
        }
        with pytest.raises(KeyError):
            TargetBody.from_isd(isd)

    def test_inconsistent_radii_raises(self):
        """ISD ``radii`` disagrees with ``BODY<code>_RADII`` by > 1 m —
        indicates stale SPICE blob or corrupted cube. Refuse to guess."""
        isd = {
            "radii": {"semimajor": 3396.19, "semiminor": 3376.2, "unit": "km"},
            "naif_keywords": {
                "BODY_CODE": 499,
                # Off by 10 km — way above the 1 m tolerance.
                "BODY499_RADII": [3400.0, 3400.0, 3376.2],
            },
        }
        with pytest.raises(ValueError, match="disagree"):
            TargetBody.from_isd(isd)


# ------------------------------------------------------------------
# mapping_to_crs: body-agnostic hardening


class TestMappingToCrsNoMarsDefault:
    def test_mars_mapping_still_works(self):
        """The common Mars path must still build a CRS when EquatorialRadius
        is explicit. Same as 0.7.x behavior for Mars files."""
        mapping = {
            "ProjectionName": "Equirectangular",
            "EquatorialRadius": 3396190.0,
            "PolarRadius": 3376200.0,
            "CenterLatitude": 0.0,
            "CenterLongitude": 0.0,
        }
        crs = mapping_to_crs(mapping)
        assert crs is not None
        # CRS ellipsoid should carry Mars radii
        ellipsoid = crs.ellipsoid
        assert abs(ellipsoid.semi_major_metre - 3396190.0) < 1.0

    def test_europa_mapping_works(self):
        """A Europa MAP file must produce a Europa CRS, not a Mars CRS."""
        mapping = {
            "ProjectionName": "Equirectangular",
            "TargetName": "Europa",
            "EquatorialRadius": 1562600.0,
            "PolarRadius": 1559500.0,
            "CenterLatitude": 0.0,
            "CenterLongitude": 0.0,
        }
        crs = mapping_to_crs(mapping)
        # Must carry Europa's radii, not Mars
        ellipsoid = crs.ellipsoid
        assert abs(ellipsoid.semi_major_metre - 1562600.0) < 1.0
        assert abs(ellipsoid.semi_major_metre - 3396190.0) > 1000

    def test_missing_equatorial_radius_raises(self):
        """Regression test for the body-agnostic refactor: a Mapping group
        with no EquatorialRadius previously silently got Mars radii. Now
        it must raise a clear error.

        Before 0.8.0 this exact input produced a Mars CRS; callers had no
        way to detect the silent fallback. The fix makes the failure
        mode loud."""
        mapping = {
            "ProjectionName": "Equirectangular",
            "TargetName": "Europa",
            "CenterLatitude": 0.0,
            "CenterLongitude": 0.0,
        }
        with pytest.raises(ValueError, match="EquatorialRadius"):
            mapping_to_crs(mapping)
