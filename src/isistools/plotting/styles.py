"""Visual styles for control network visualization.

Fixes the Qmos/Qnet problem where registered points are nearly
invisible while unregistered points are prominent. Here, registered
points are bright and clear, unregistered are secondary, and ignored
points are dimmed.
"""

from __future__ import annotations

# Point styles keyed by status (as classified by io.controlnet)
CNET_POINT_STYLES = {
    "registered": {
        "color": "#2ecc71",  # bright green — the main thing you want to see
        "size": 8,
        "alpha": 0.9,
        "marker": "circle",
        "line_color": "#27ae60",
        "line_width": 1,
    },
    "unregistered": {
        "color": "#e74c3c",  # red — needs attention but secondary
        "size": 6,
        "alpha": 0.7,
        "marker": "triangle",
        "line_color": "#c0392b",
        "line_width": 1,
    },
    "ignored": {
        "color": "#95a5a6",  # gray — deliberately faded
        "size": 4,
        "alpha": 0.4,
        "marker": "x",
        "line_color": "#7f8c8d",
        "line_width": 1,
    },
    "selected": {
        "color": "#f1c40f",  # yellow — for interactive selection highlights
        "size": 12,
        "alpha": 1.0,
        "marker": "star",
        "line_color": "#f39c12",
        "line_width": 2,
    },
}

# Color map for status values (used in hvplot color mapping)
STATUS_COLOR_MAP = {
    "registered": "#2ecc71",
    "unregistered": "#e74c3c",
    "ignored": "#95a5a6",
}

# Footprint styles
FOOTPRINT_STYLES = {
    "default": {
        "fill_color": "#3498db",
        "fill_alpha": 0.15,
        "line_color": "#2980b9",
        "line_width": 1.5,
    },
    "selected": {
        "fill_color": "#e74c3c",
        "fill_alpha": 0.3,
        "line_color": "#c0392b",
        "line_width": 2.5,
    },
    "hover": {
        "fill_color": "#f39c12",
        "fill_alpha": 0.25,
        "line_color": "#e67e22",
        "line_width": 2.0,
    },
}

# Image display defaults
IMAGE_DEFAULTS = {
    "cmap": "gray",
    "clim_percentile": (1, 99),  # stretch to 1st–99th percentile
}


def status_to_bokeh_style(status: str) -> dict:
    """Get Bokeh-compatible glyph style for a point status.

    Returns a dict suitable for passing to Bokeh scatter kwargs.
    """
    style = CNET_POINT_STYLES.get(status, CNET_POINT_STYLES["unregistered"])
    return {
        "fill_color": style["color"],
        "fill_alpha": style["alpha"],
        "size": style["size"],
        "marker": style["marker"],
        "line_color": style["line_color"],
        "line_width": style["line_width"],
    }
