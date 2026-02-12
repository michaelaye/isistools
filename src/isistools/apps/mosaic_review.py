"""Mosaic Review — Qmos replacement.

Interactive Panel app combining:
- Footprint overview map (click to select)
- Rasterized image viewer for selected cube
- Control network overlay with proper styling
- Point info panel

Works in notebooks (via ``.panel()``) and as a standalone server
(via ``.serve()`` or the CLI).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import holoviews as hv
import panel as pn

from isistools.apps.components import CnetInfoPanel, CnetSelector, CubeListSelector
from isistools.io.controlnet import load_cnet
from isistools.io.cubes import load_cube, read_label
from isistools.io.footprints import load_footprints, read_cube_list
from isistools.plotting.cnet_overlay import cnet_to_geodataframe
from isistools.plotting.footprint_map import footprint_map, footprint_map_with_cnet
from isistools.plotting.image_viewer import image_with_cnet, image_plot

if TYPE_CHECKING:
    import geopandas as gpd
    import pandas as pd

pn.extension("tabulator")
hv.extension("bokeh")


class MosaicReview:
    """Interactive mosaic review application.

    Combines a footprint overview map with an image browser and
    control network overlay. Replaces Qmos.

    Parameters
    ----------
    cube_list : path-like, optional
        Path to a cube list file, or list of cube paths.
    cnet_path : path-like, optional
        Path to a control network file.

    Examples
    --------
    In a notebook::

        app = MosaicReview("cubes.lis", cnet_path="control.net")
        app.panel()

    From CLI (opens browser)::

        app = MosaicReview("cubes.lis")
        app.serve()
    """

    def __init__(
        self,
        cube_list: str | Path | list[str | Path] | None = None,
        cnet_path: str | Path | None = None,
    ):
        self._cube_paths: list[Path] = []
        self._footprints: gpd.GeoDataFrame | None = None
        self._cnet_df: pd.DataFrame | None = None
        self._cnet_gdf: gpd.GeoDataFrame | None = None
        self._selected_cube: Path | None = None

        # Widgets
        self._cube_selector = CubeListSelector(default_path=cube_list)
        self._cnet_selector = CnetSelector(default_path=cnet_path)
        self._cnet_info = CnetInfoPanel()

        # Plot panes
        self._map_pane = pn.pane.HoloViews(hv.Div("Load a cube list to begin"), sizing_mode="stretch_both")
        self._image_pane = pn.pane.HoloViews(hv.Div("Select an image from the map"), sizing_mode="stretch_both")

        # Image selector dropdown (populated after loading cubes)
        self._image_dropdown = pn.widgets.Select(
            name="Image", options=[], width=400,
        )
        self._image_dropdown.param.watch(self._on_image_selected, "value")

        # Wire up callbacks
        self._cube_selector.on_load(self._on_cubes_loaded)
        self._cnet_selector.on_load(self._on_cnet_loaded)

        # Auto-load if paths provided
        if cube_list is not None:
            self._auto_load_cubes(cube_list)
        if cnet_path is not None:
            self._auto_load_cnet(cnet_path)

    def _auto_load_cubes(self, cube_list):
        """Load cubes at init time (non-interactive)."""
        try:
            if isinstance(cube_list, (str, Path)):
                path = Path(cube_list)
                if path.is_file():
                    self._cube_paths = read_cube_list(path)
                else:
                    self._cube_paths = [path]
            else:
                self._cube_paths = [Path(p) for p in cube_list]

            self._footprints = load_footprints(self._cube_paths, skip_errors=True)
            self._update_map()
            self._update_image_dropdown()
        except Exception as e:
            self._map_pane.object = hv.Div(f"<b>Error loading cubes:</b> {e}")

    def _auto_load_cnet(self, cnet_path):
        """Load control network at init time."""
        try:
            self._cnet_df = load_cnet(cnet_path)
            self._cnet_gdf = cnet_to_geodataframe(self._cnet_df)
            self._cnet_info.update(self._cnet_df)
            self._update_map()
        except Exception as e:
            self._cnet_info._content.object = f"*Error: {e}*"

    def _on_cubes_loaded(self, cube_paths: list[Path]):
        """Callback when cube list is loaded via widget."""
        self._cube_paths = cube_paths
        self._footprints = load_footprints(cube_paths, skip_errors=True)
        self._update_map()
        self._update_image_dropdown()

    def _on_cnet_loaded(self, cnet_df):
        """Callback when control network is loaded via widget."""
        self._cnet_df = cnet_df
        self._cnet_gdf = cnet_to_geodataframe(cnet_df)
        self._cnet_info.update(cnet_df)
        self._update_map()

    def _update_map(self):
        """Refresh the footprint map."""
        if self._footprints is None or self._footprints.empty:
            self._map_pane.object = hv.Div("No footprints loaded")
            return

        if self._cnet_gdf is not None and not self._cnet_gdf.empty:
            self._map_pane.object = footprint_map_with_cnet(
                self._footprints, self._cnet_gdf, title="Mosaic Footprints"
            )
        else:
            self._map_pane.object = footprint_map(
                self._footprints, title="Mosaic Footprints"
            )

    def _update_image_dropdown(self):
        """Populate the image selector with loaded cubes."""
        options = {p.name: str(p) for p in self._cube_paths}
        self._image_dropdown.options = options

    def _on_image_selected(self, event):
        """Show the selected image."""
        if not event.new:
            return
        cube_path = Path(event.new)
        try:
            da = load_cube(cube_path)
            if self._cnet_df is not None:
                # Match via spacecraft clock count
                clock_lookup = _build_serial_lookup([cube_path])
                label = read_label(cube_path)
                inst = label["IsisCube"]["Instrument"]
                clock = str(
                    inst.get("SpacecraftClockCount",
                             inst.get("SpacecraftClockStartCount", ""))
                )
                # Find the serial number that contains this clock count
                matching = [
                    sn for sn in self._cnet_df["serialnumber"].unique()
                    if clock and clock in sn
                ]
                if matching:
                    sn = matching[0]
                    self._image_pane.object = image_with_cnet(da, self._cnet_df, serial_number=sn)
                else:
                    self._image_pane.object = image_plot(da)
            else:
                self._image_pane.object = image_plot(da)
        except Exception as e:
            self._image_pane.object = hv.Div(f"<b>Error:</b> {e}")

    def panel(self) -> pn.viewable.Viewable:
        """Return the Panel layout for notebook or server use.

        Returns
        -------
        panel.viewable.Viewable
            The complete app layout.
        """
        header = pn.Row(
            self._cube_selector,
            self._cnet_selector,
        )
        controls = pn.Row(self._image_dropdown)

        main = pn.Row(
            self._map_pane,
            pn.Column(controls, self._image_pane),
            sizing_mode="stretch_both",
        )

        sidebar = pn.Column(
            self._cnet_info,
            width=300,
        )

        return pn.Column(
            pn.pane.Markdown("# Mosaic Review"),
            header,
            pn.Row(main, sidebar, sizing_mode="stretch_both"),
            sizing_mode="stretch_both",
        )

    def serve(self, port: int = 0, show: bool = True, **kwargs):
        """Launch as a standalone Panel server in the browser.

        Parameters
        ----------
        port : int
            Port number. 0 means auto-select.
        show : bool
            Open browser automatically.
        **kwargs
            Additional arguments passed to ``panel.serve``.
        """
        pn.serve(self.panel(), port=port, show=show, title="Mosaic Review", **kwargs)
