"""Tests for AlphaCube handling in csm2map.

Regression guard for a 2026-09-30 defect: csm2map projected every cropped cube
as though it began at line 1 of its parent, because ALE builds the ISD from the
cube's Core dimensions and SPICE (both describing the parent) while csm2map
applied the resulting model to the crop's own pixel grid.

Measured on ``G01_018621_1678_XI_12S339W`` before the fix: the same cube cropped
at lines 2000, 8000 and 20000 projected to an identical scene centre, putting it
11.9 km, 47.7 km and 119.2 km from where ``campt`` placed it. After the fix, a
crop at line 17000 of ``F20_043544_1795_XI_00S359W`` agreed with ``cam2map`` to
0.6 pixels (3 m at 5 m/px), measured by phase correlation.
"""

import pytest

from isistools.csm2map.camera import AlphaCube, CroppedRasterGM, read_alpha_cube


class _Coord:
    """Stand-in for csmapi.ImageCoord."""

    def __init__(self, line, samp):
        self.line = line
        self.samp = samp


class _StubModel:
    """Minimal RasterGM stand-in recording what it was asked."""

    def __init__(self):
        self.last_image_pt = None

    def groundToImage(self, _ground_pt):  # noqa: N802 - mirrors csmapi
        return _Coord(17299.0, 2500.0)

    def imageToGround(self, image_pt, _height=0.0):  # noqa: N802
        self.last_image_pt = (image_pt.line, image_pt.samp)
        return "ground"

    def somethingElse(self):  # noqa: N802
        return "delegated"


class TestReadAlphaCube:
    def test_absent_group_returns_none(self, monkeypatch):
        monkeypatch.setattr(
            "isistools.io.cubes.read_label",
            lambda _p: {"IsisCube": {"Core": {"Dimensions": {"Lines": 10, "Samples": 5}}}},
        )
        assert read_alpha_cube("x.cub") is None

    def test_identity_crop_returns_none(self, monkeypatch):
        """A cube cropped from its own first pixel needs no correction."""
        monkeypatch.setattr(
            "isistools.io.cubes.read_label",
            lambda _p: {
                "IsisCube": {
                    "Core": {"Dimensions": {"Lines": 600, "Samples": 5000}},
                    "AlphaCube": {"AlphaStartingLine": 0.5, "AlphaStartingSample": 0.5},
                }
            },
        )
        assert read_alpha_cube("x.cub") is None

    def test_offsets_use_isis_edge_convention(self, monkeypatch):
        """ISIS records leading edges; CSM puts the first pixel centre at 0.5.

        ``crop line=17000`` writes ``AlphaStartingLine = 16999.5``, and the
        crop's line 1 is the parent's line 17000, so the offset is 16999.
        """
        monkeypatch.setattr(
            "isistools.io.cubes.read_label",
            lambda _p: {
                "IsisCube": {
                    "Core": {"Dimensions": {"Lines": 600, "Samples": 5000}},
                    "AlphaCube": {
                        "AlphaStartingLine": 16999.5,
                        "AlphaStartingSample": 0.5,
                    },
                }
            },
        )
        ac = read_alpha_cube("x.cub")
        assert ac == AlphaCube(
            line_offset=16999.0, sample_offset=0.0, lines=600, samples=5000
        )

    def test_malformed_group_returns_none(self, monkeypatch):
        monkeypatch.setattr(
            "isistools.io.cubes.read_label",
            lambda _p: {
                "IsisCube": {
                    "Core": {"Dimensions": {"Lines": 600, "Samples": 5000}},
                    "AlphaCube": {"AlphaStartingLine": "not-a-number"},
                }
            },
        )
        assert read_alpha_cube("x.cub") is None


class TestCroppedRasterGM:
    @pytest.fixture
    def wrapped(self):
        alpha = AlphaCube(
            line_offset=16999.0, sample_offset=0.0, lines=600, samples=5000
        )
        stub = _StubModel()
        return CroppedRasterGM(stub, alpha), stub

    def test_image_to_ground_shifts_into_parent_frame(self, wrapped):
        model, stub = wrapped
        model.imageToGround(_Coord(300.0, 2500.0), 0.0)
        assert stub.last_image_pt == (17299.0, 2500.0)

    def test_ground_to_image_shifts_back_into_crop_frame(self, wrapped, monkeypatch):
        model, _ = wrapped
        monkeypatch.setitem(__import__("sys").modules, "csmapi", _FakeCsmapi())
        ip = model.groundToImage(object())
        assert (ip.line, ip.samp) == (300.0, 2500.0)

    def test_round_trip_is_identity(self, wrapped, monkeypatch):
        model, _ = wrapped
        monkeypatch.setitem(__import__("sys").modules, "csmapi", _FakeCsmapi())
        ip = model.groundToImage(object())
        model.imageToGround(ip, 0.0)
        assert model._model.last_image_pt == (17299.0, 2500.0)

    def test_unknown_attributes_delegate(self, wrapped):
        model, _ = wrapped
        assert model.somethingElse() == "delegated"


class _FakeCsmapi:
    ImageCoord = _Coord

    class ImageVector:
        def __init__(self, line, samp):
            self.line = line
            self.samp = samp
