"""Tiepoint Review — Qnet replacement.

Interactive Panel app for reviewing control network tie points
between image pairs:
- Side-by-side image pair viewer with correct (north-up) orientation
- Control network points with improved color scheme
- Residual vector display
- Point detail panel

Key improvements over Qnet:
- Images shown in map orientation, not detector readout order
- Linked crosshairs between image pair
- Proper point coloring (registered=visible, unregistered=secondary)
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import holoviews as hv
import panel as pn

from isistools.apps.components import CnetInfoPanel, PointDetailPanel
from isistools.io.controlnet import load_cnet
from isistools.io.cubes import load_cube, read_label
from isistools.io.footprints import read_cube_list
from isistools.plotting.cnet_overlay import cnet_points_image, cnet_residual_vectors
from isistools.plotting.image_viewer import image_plot

if TYPE_CHECKING:
    import pandas as pd

pn.extension("tabulator")
hv.extension("bokeh")


def _find_image_pairs(cnet_df: pd.DataFrame) -> list[tuple[str, str]]:
    """Find all unique image pairs that share control points.

    Returns a list of (serial1, serial2) tuples, sorted by the number
    of shared points (most shared first).
    """
    import pandas as pd

    # For each point, get all serial numbers
    point_serials = cnet_df.groupby("pointId")["serialnumber"].apply(set)

    pair_counts: dict[tuple[str, str], int] = {}
    for serials in point_serials:
        serials_list = sorted(serials)
        for i, s1 in enumerate(serials_list):
            for s2 in serials_list[i + 1:]:
                pair = (s1, s2)
                pair_counts[pair] = pair_counts.get(pair, 0) + 1

    # Sort by count descending
    sorted_pairs = sorted(pair_counts.items(), key=lambda x: -x[1])
    return [pair for pair, count in sorted_pairs]


from isistools.io.cubes import build_serial_lookup, match_serials_to_cubes


class TiepointReview:
    """Interactive tie point review application.

    Replaces Qnet with proper image orientation and improved
    point visualization.

    Parameters
    ----------
    cube_list : path-like
        Path to cube list file or list of cube paths.
    cnet_path : path-like
        Path to control network file.
    """

    def __init__(
        self,
        cube_list: str | Path | list[str | Path],
        cnet_path: str | Path,
    ):
        # Load data
        if isinstance(cube_list, (str, Path)):
            path = Path(cube_list)
            if path.suffix in (".lis", ".txt", ".list", ""):
                self._cube_paths = read_cube_list(path)
            else:
                self._cube_paths = [path]
        else:
            self._cube_paths = [Path(p) for p in cube_list]

        self._cnet_df = load_cnet(cnet_path)
        self._cnet_info = CnetInfoPanel()
        self._cnet_info.update(self._cnet_df)
        self._point_detail = PointDetailPanel()

        # Find image pairs
        self._pairs = _find_image_pairs(self._cnet_df)
        self._serial_to_path = match_serials_to_cubes(
            self._cnet_df["serialnumber"].unique().tolist(),
            self._cube_paths,
        )

        # Build pair labels for selector
        pair_labels = {}
        for s1, s2 in self._pairs:
            p1 = self._serial_to_path.get(s1)
            p2 = self._serial_to_path.get(s2)
            name1 = p1.stem if p1 else s1.split("/")[-1]
            name2 = p2.stem if p2 else s2.split("/")[-1]
            # Count shared points
            shared = self._cnet_df[
                self._cnet_df["serialnumber"].isin([s1, s2])
            ]["pointId"].value_counts()
            n_shared = (shared >= 2).sum()
            label = f"{name1} ↔ {name2} ({n_shared} pts)"
            pair_labels[label] = (s1, s2)

        self._pair_labels = pair_labels

        # Widgets
        self._pair_selector = pn.widgets.Select(
            name="Image Pair",
            options=list(pair_labels.keys()),
            width=500,
        )
        self._show_residuals = pn.widgets.Checkbox(
            name="Show Residual Vectors", value=False,
        )
        self._residual_scale = pn.widgets.FloatSlider(
            name="Residual Scale",
            start=1, end=100, value=10, step=1,
            width=200,
        )

        # Plot panes
        self._left_pane = pn.pane.HoloViews(
            hv.Div("Select an image pair"), sizing_mode="stretch_both"
        )
        self._right_pane = pn.pane.HoloViews(
            hv.Div(""), sizing_mode="stretch_both"
        )

        # Wire up
        self._pair_selector.param.watch(self._on_pair_selected, "value")
        self._show_residuals.param.watch(self._on_pair_selected, "value")
        self._residual_scale.param.watch(self._on_pair_selected, "value")

        # Auto-select first pair
        if self._pair_labels:
            self._on_pair_selected(None)

    def _on_pair_selected(self, event):
        """Load and display the selected image pair."""
        label = self._pair_selector.value
        if label not in self._pair_labels:
            return

        s1, s2 = self._pair_labels[label]
        path1 = self._serial_to_path.get(s1)
        path2 = self._serial_to_path.get(s2)

        if path1 is None or path2 is None:
            msg = "Could not find cube files for this pair"
            self._left_pane.object = hv.Div(f"<b>{msg}</b>")
            self._right_pane.object = hv.Div("")
            return

        try:
            da1 = load_cube(path1)
            da2 = load_cube(path2)

            # Build left image + cnet
            img1 = image_plot(da1, title=path1.stem)
            pts1 = cnet_points_image(self._cnet_df, serial_number=s1)
            left = img1 * pts1

            if self._show_residuals.value:
                vecs1 = cnet_residual_vectors(
                    self._cnet_df, serial_number=s1,
                    scale=self._residual_scale.value,
                )
                left = left * vecs1

            # Build right image + cnet
            img2 = image_plot(da2, title=path2.stem)
            pts2 = cnet_points_image(self._cnet_df, serial_number=s2)
            right = img2 * pts2

            if self._show_residuals.value:
                vecs2 = cnet_residual_vectors(
                    self._cnet_df, serial_number=s2,
                    scale=self._residual_scale.value,
                )
                right = right * vecs2

            self._left_pane.object = left
            self._right_pane.object = right

        except Exception as e:
            self._left_pane.object = hv.Div(f"<b>Error:</b> {e}")
            self._right_pane.object = hv.Div("")

    def panel(self) -> pn.viewable.Viewable:
        """Return the Panel layout.

        Returns
        -------
        panel.viewable.Viewable
        """
        controls = pn.Row(
            self._pair_selector,
            self._show_residuals,
            self._residual_scale,
        )

        images = pn.Row(
            self._left_pane,
            self._right_pane,
            sizing_mode="stretch_both",
        )

        sidebar = pn.Column(
            self._cnet_info,
            pn.layout.Divider(),
            self._point_detail,
            width=300,
        )

        return pn.Column(
            pn.pane.Markdown("# Tiepoint Review"),
            controls,
            pn.Row(images, sidebar, sizing_mode="stretch_both"),
            sizing_mode="stretch_both",
        )

    def serve(self, port: int = 0, show: bool = True, **kwargs):
        """Launch as a standalone Panel server.

        Parameters
        ----------
        port : int
            Port number (0 = auto).
        show : bool
            Open browser.
        """
        pn.serve(
            self.panel(), port=port, show=show,
            title="Tiepoint Review", **kwargs,
        )
