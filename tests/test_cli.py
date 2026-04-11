"""Tests for the Typer CLI wiring in isistools.cli.

These tests focus on command *registration* and *import correctness* — they
verify that each command can be invoked without blowing up on import-time or
body-level NameError, without requiring the heavy external deps (ISIS,
geopandas data files, etc.) to actually be present on the system under test.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from isistools.cli import app

runner = CliRunner()


def test_cli_module_imports():
    """Importing isistools.cli must not raise.

    This is the broadest possible smoke check — it catches syntax errors,
    unresolved top-level imports, and any other module-load-time crash that
    would take the entire CLI down.
    """
    import isistools.cli  # noqa: F401


def test_overlaps_help():
    """`isistools overlaps --help` must exit 0.

    Typer introspects the function signature to render help, so this fails
    only if the command fails to register or the option types don't resolve.
    """
    result = runner.invoke(app, ["overlaps", "--help"])
    assert result.exit_code == 0, result.output
    assert "findimageoverlaps" in result.output


def test_overlaps_png_path_does_not_hit_undefined_gpd():
    """Regression test for F821: `overlaps --png` referenced `gpd` without import.

    The ``overlaps`` command's PNG-plot branch calls
    ``gpd.GeoDataFrame(...).plot(...)`` inside its per-row loop (cli.py line
    ~395). Before the fix (commit before 0.7.0), ``gpd`` was never imported —
    so any user who ran ``isistools overlaps <cubelist> --png`` hit
    ``NameError: name 'gpd' is not defined`` the moment the first row was
    drawn. The bug survived 0.6.0 because ``ruff`` wasn't gating merges and
    no test exercised the PNG path.

    This test mocks out the ISIS ``findimageoverlaps`` subprocess and the
    ``parse_overlap_list`` helper so the PNG branch can execute without any
    external tooling, then asserts the command returns exit-code 0. Before
    the fix, the assertion fails with a NameError in ``result.exception``.
    """
    # Import here so the test fails clearly if geopandas isn't installed,
    # rather than at collection time.
    import geopandas as gpd  # type: ignore[import-not-found]
    from shapely.geometry import Polygon

    # Build a minimal GeoDataFrame mimicking parse_overlap_list's output.
    fake_gdf = gpd.GeoDataFrame(
        {
            "serials": ["MRO/CTX/A,MRO/CTX/B"],
            "zone_type": ["2-way overlap"],
            "area_deg2": [0.123],
            "n_images": [2],
            "images": [["a.cub", "b.cub"]],
            "geometry": [Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        },
        crs="EPSG:4326",
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        cubelist = Path(tmpdir) / "cubes.lis"
        cubelist.write_text("/fake/cube.cub\n")
        png_out = Path(tmpdir) / "overlaps.png"

        # Force a non-interactive matplotlib backend for CI/headless runs.
        import matplotlib

        matplotlib.use("Agg", force=True)

        fake_proc = MagicMock(returncode=0, stderr="", stdout="")
        with (
            patch("subprocess.run", return_value=fake_proc),
            patch(
                "isistools.io.overlaps.parse_overlap_list",
                return_value=fake_gdf,
            ),
        ):
            result = runner.invoke(
                app,
                [
                    "overlaps",
                    str(cubelist),
                    "--png",
                    "--png-path",
                    str(png_out),
                ],
            )

        if result.exit_code != 0:
            # Surface the NameError (or whatever else) so the failure is legible.
            raise AssertionError(
                f"overlaps --png failed with exit_code={result.exit_code}\n"
                f"stdout:\n{result.output}\n"
                f"exception: {result.exception!r}"
            ) from result.exception

        assert png_out.exists(), "PNG output was not produced"


@pytest.mark.parametrize(
    "command",
    [
        "mosaic",
        "tiepoints",
        "footprints",
        "footprintinit",
        "spiceinit",
        "overlaps",
        "cnet-info",
        "csm2map",
        "csm2map-compare",
    ],
)
def test_every_command_has_help(command):
    """Every registered command must render `--help` without crashing."""
    result = runner.invoke(app, [command, "--help"])
    assert result.exit_code == 0, f"{command} --help failed: {result.output}"
