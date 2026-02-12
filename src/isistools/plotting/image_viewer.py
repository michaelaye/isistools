"""Rasterized image viewer for ISIS cubes.

Uses hvplot with datashader rasterization for responsive display
of large planetary images (e.g., full HiRISE strips). Supports
optional control network point overlays.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import holoviews as hv
import hvplot.xarray  # noqa: F401
import numpy as np

from isistools.plotting.styles import IMAGE_DEFAULTS

if TYPE_CHECKING:
    import pandas as pd
    import xarray as xr

hv.extension("bokeh")


def image_plot(
    da: xr.DataArray,
    rasterize: bool = True,
    cmap: str | None = None,
    percentile_stretch: tuple[float, float] | None = None,
    title: str | None = None,
    width: int = 600,
    height: int = 600,
    responsive: bool = False,
) -> hv.Element:
    """Display an ISIS cube image with rasterized rendering.

    Parameters
    ----------
    da : xr.DataArray
        Image data (as returned by :func:`isistools.io.cubes.load_cube`).
    rasterize : bool
        If True, use datashader for rasterization. Essential for
        large images (HiRISE, CTX).
    cmap : str, optional
        Colormap name. Defaults to 'gray'.
    percentile_stretch : tuple of float, optional
        Percentile range for contrast stretch, e.g. (1, 99).
        Defaults to styles.IMAGE_DEFAULTS["clim_percentile"].
    title : str, optional
        Plot title. Defaults to the cube filename.
    width, height : int
        Plot dimensions.
    responsive : bool
        If True, plot fills available space.

    Returns
    -------
    holoviews.Element
        Interactive image plot.
    """
    if cmap is None:
        cmap = IMAGE_DEFAULTS["cmap"]
    if percentile_stretch is None:
        percentile_stretch = IMAGE_DEFAULTS["clim_percentile"]

    if title is None:
        title = da.attrs.get("cube_path", "Image")
        # Just show filename
        if "/" in title:
            title = title.rsplit("/", 1)[-1]

    # Compute contrast limits from a subsample for performance
    clim = _compute_clim(da, percentile_stretch)

    def _deduplicate_tools(plot, element):
        """Remove duplicate toolbar buttons added by Bokeh defaults."""
        seen = set()
        unique = []
        for tool in plot.state.toolbar.tools:
            name = type(tool).__name__
            if name not in seen:
                seen.add(name)
                unique.append(tool)
        plot.state.toolbar.tools = unique

    plot_kwargs = dict(
        rasterize=rasterize,
        cmap=cmap,
        clim=clim,
        title=title,
        tools=["wheel_zoom", "pan", "reset", "crosshair"],
        aspect="equal",
    )
    if responsive:
        plot_kwargs["responsive"] = True
    else:
        plot_kwargs["width"] = width
        plot_kwargs["height"] = height

    return da.hvplot.image(**plot_kwargs).opts(hooks=[_deduplicate_tools])


def image_pair_plot(
    da_left: xr.DataArray,
    da_right: xr.DataArray,
    link_axes: bool = True,
    **kwargs,
) -> hv.Layout:
    """Side-by-side image pair viewer with optional linked axes.

    Designed for comparing overlapping images in a control network,
    similar to Qnet's pair view but with correct (north-up) orientation.

    Parameters
    ----------
    da_left, da_right : xr.DataArray
        The two images to compare.
    link_axes : bool
        If True, pan/zoom is synchronized between the two images.
    **kwargs
        Additional arguments passed to :func:`image_plot`.

    Returns
    -------
    holoviews.Layout
        Side-by-side plot layout.
    """
    left = image_plot(da_left, **kwargs)
    right = image_plot(da_right, **kwargs)

    if link_axes:
        # Link the x/y ranges so pan/zoom is synchronized
        left = left.opts(shared_axes=True)
        right = right.opts(shared_axes=True)

    return left + right


def image_with_cnet(
    da: xr.DataArray,
    cnet_df: pd.DataFrame,
    serial_number: str | None = None,
    **kwargs,
) -> hv.Element:
    """Image plot with control network points overlaid.

    Parameters
    ----------
    da : xr.DataArray
        Image data.
    cnet_df : pd.DataFrame
        Control network (as returned by :func:`isistools.io.controlnet.load_cnet`).
    serial_number : str, optional
        If provided, only show measures matching this serial number.
        Otherwise shows all measures.
    **kwargs
        Passed to :func:`image_plot`.

    Returns
    -------
    holoviews.Element
        Image with cnet points overlay.
    """
    from isistools.plotting.cnet_overlay import cnet_points_image

    img = image_plot(da, **kwargs)
    points = cnet_points_image(cnet_df, serial_number=serial_number)

    return img * points


def _compute_clim(
    da: xr.DataArray,
    percentile: tuple[float, float],
    max_samples: int = 500_000,
) -> tuple[float, float]:
    """Compute contrast stretch limits from a data subsample.

    Subsamples the data for performance on large arrays.
    Ignores NaN/special pixel values.
    """
    data = da.values.ravel()
    # Remove NaN and ISIS special pixels (very large negative values)
    valid = data[np.isfinite(data) & (data > -1e30)]

    if len(valid) == 0:
        return (0, 1)

    if len(valid) > max_samples:
        rng = np.random.default_rng(42)
        valid = rng.choice(valid, max_samples, replace=False)

    lo = np.percentile(valid, percentile[0])
    hi = np.percentile(valid, percentile[1])

    return (float(lo), float(hi))
