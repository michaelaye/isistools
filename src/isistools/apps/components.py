"""Reusable Panel widget components for isistools apps.

Provides common UI elements: file selectors, info panels,
status displays, etc.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import panel as pn

if TYPE_CHECKING:
    import geopandas as gpd
    import pandas as pd

pn.extension("tabulator")


class CubeListSelector(pn.viewable.Viewer):
    """Widget for selecting a cube list file and loading it."""

    def __init__(self, default_path: str | Path | None = None, **params):
        super().__init__(**params)
        self._path_input = pn.widgets.TextInput(
            name="Cube List",
            value=str(default_path) if default_path else "",
            placeholder="/path/to/cubelist.lis",
            width=400,
        )
        self._load_btn = pn.widgets.Button(
            name="Load", button_type="primary", width=80
        )
        self._status = pn.pane.Str("", styles={"color": "#666"})

        self._on_load_callbacks: list = []
        self._load_btn.on_click(self._on_load)

    def on_load(self, callback):
        """Register a callback for when a cube list is loaded.

        Callback receives a list of Path objects.
        """
        self._on_load_callbacks.append(callback)

    def _on_load(self, event):
        from isistools.io.footprints import read_cube_list

        path = Path(self._path_input.value)
        if not path.exists():
            self._status.object = f"File not found: {path}"
            return
        try:
            cubes = read_cube_list(path)
            self._status.object = f"Loaded {len(cubes)} cubes"
            for cb in self._on_load_callbacks:
                cb(cubes)
        except Exception as e:
            self._status.object = f"Error: {e}"

    def __panel__(self):
        return pn.Row(self._path_input, self._load_btn, self._status)


class CnetSelector(pn.viewable.Viewer):
    """Widget for selecting and loading a control network file."""

    def __init__(self, default_path: str | Path | None = None, **params):
        super().__init__(**params)
        self._path_input = pn.widgets.TextInput(
            name="Control Net",
            value=str(default_path) if default_path else "",
            placeholder="/path/to/control.net",
            width=400,
        )
        self._load_btn = pn.widgets.Button(
            name="Load", button_type="primary", width=80
        )
        self._status = pn.pane.Str("", styles={"color": "#666"})

        self._on_load_callbacks: list = []
        self._load_btn.on_click(self._on_load)

    def on_load(self, callback):
        """Register a callback for when a cnet is loaded.

        Callback receives a pd.DataFrame.
        """
        self._on_load_callbacks.append(callback)

    def _on_load(self, event):
        from isistools.io.controlnet import load_cnet

        path = Path(self._path_input.value)
        if not path.exists():
            self._status.object = f"File not found: {path}"
            return
        try:
            df = load_cnet(path)
            self._status.object = f"Loaded {len(df)} measures"
            for cb in self._on_load_callbacks:
                cb(df)
        except Exception as e:
            self._status.object = f"Error: {e}"

    def __panel__(self):
        return pn.Row(self._path_input, self._load_btn, self._status)


class CnetInfoPanel(pn.viewable.Viewer):
    """Displays summary statistics for a control network."""

    def __init__(self, **params):
        super().__init__(**params)
        self._content = pn.pane.Markdown("*No control network loaded*")

    def update(self, cnet_df: pd.DataFrame):
        from isistools.io.controlnet import cnet_summary

        stats = cnet_summary(cnet_df)
        md = f"""### Control Network Summary
| Metric | Value |
|--------|-------|
| Points | {stats['n_points']} |
| Measures | {stats['n_measures']} |
| Images | {stats['n_images']} |
| Registered | {stats['n_registered']} |
| Unregistered | {stats['n_unregistered']} |
| Ignored | {stats['n_ignored']} |
| Mean Residual | {stats['mean_residual']:.4f} |
| Max Residual | {stats['max_residual']:.4f} |
"""
        self._content.object = md

    def __panel__(self):
        return self._content


class PointDetailPanel(pn.viewable.Viewer):
    """Shows details for a selected control point."""

    def __init__(self, **params):
        super().__init__(**params)
        self._content = pn.pane.Markdown("*Click a point to see details*")

    def update(self, point_id: str, cnet_df: pd.DataFrame):
        measures = cnet_df[cnet_df["pointId"] == point_id]
        if measures.empty:
            self._content.object = f"*Point {point_id} not found*"
            return

        first = measures.iloc[0]
        lines = [f"### Point: {point_id}"]
        lines.append(f"**Type:** {first.get('pointType', 'N/A')}")
        lines.append(f"**Measures:** {len(measures)}")
        lines.append("")
        lines.append("| Image | Sample | Line | Residual | Status |")
        lines.append("|-------|--------|------|----------|--------|")

        for _, m in measures.iterrows():
            sn = m.get("serialnumber", "?")
            # Show just the last part of the serial number for readability
            sn_short = sn.rsplit("/", 1)[-1] if "/" in str(sn) else str(sn)
            lines.append(
                f"| {sn_short} | {m.get('sample', 0):.1f} | "
                f"{m.get('line', 0):.1f} | "
                f"{m.get('residual_magnitude', 0):.3f} | "
                f"{m.get('status', '?')} |"
            )

        self._content.object = "\n".join(lines)

    def __panel__(self):
        return self._content
