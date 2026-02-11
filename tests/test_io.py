"""Tests for isistools.io modules."""

from isistools.io.controlnet import _classify_point_status


class TestClassifyPointStatus:
    def test_ignored_point(self):
        row = {"pointIgnore": True, "measureIgnore": False, "measureType": 0}
        import pandas as pd
        assert _classify_point_status(pd.Series(row)) == "ignored"

    def test_ignored_measure(self):
        row = {"pointIgnore": False, "measureIgnore": True, "measureType": 0}
        import pandas as pd
        assert _classify_point_status(pd.Series(row)) == "ignored"

    def test_registered_by_type(self):
        row = {"pointIgnore": False, "measureIgnore": False, "measureType": 2}
        import pandas as pd
        assert _classify_point_status(pd.Series(row)) == "registered"

    def test_registered_by_residual(self):
        row = {
            "pointIgnore": False, "measureIgnore": False, "measureType": 0,
            "residualSample": 0.5, "residualLine": -0.3,
        }
        import pandas as pd
        assert _classify_point_status(pd.Series(row)) == "registered"

    def test_unregistered(self):
        row = {
            "pointIgnore": False, "measureIgnore": False, "measureType": 0,
            "residualSample": 0.0, "residualLine": 0.0,
        }
        import pandas as pd
        assert _classify_point_status(pd.Series(row)) == "unregistered"
