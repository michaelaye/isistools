---
format:
  html:
    theme: solar
---

# isistools

Python-based review tools for ISIS3 coregistration workflows. Replaces Qmos and Qnet with modern, interactive visualization using HoloViews, Panel, and datashader.

## Why?

ISIS's Qmos and Qnet have several pain points:

- **Qmos** requires level-1 unprojected images but displays footprints in map projection (confusing), and renders image content very slowly.
- **Qnet** requires level-2 map-projected images but displays them flipped in detector readout order (makes comparison difficult).
- **Both** use a color scheme where registered control points are nearly invisible while unregistered points are prominent — the opposite of what you want during review.

**isistools** fixes all of these by building on the Python geospatial stack.

## Features

- **Footprint map**: Interactive map of image footprints using geopandas + hvplot, with hover info and click-to-select.
- **Image viewer**: Fast rasterized image display via rioxarray + datashader. Images always shown in correct (north-up) orientation.
- **Control network overlay**: Tie points with sensible colors (registered = bright green, unregistered = red, ignored = gray). Residual vectors. Point detail on click.
- **Dual interface**: Same code works in Jupyter notebooks and as standalone Panel apps in the browser.
- **CLI**: `isistools mosaic`, `isistools tiepoints`, `isistools footprints` commands that launch Panel apps.

## Installation

```bash
pip install -e .
```

## Usage

### CLI

```bash
# Mosaic review (Qmos replacement)
isistools mosaic cubes.lis --cnet control.net

# Tiepoint review (Qnet replacement)
isistools tiepoints cubes.lis control.net

# Quick footprint map
isistools footprints cubes.lis

# Control network summary stats
isistools cnet-info control.net
```

### Notebook

```python
from isistools.apps.mosaic_review import MosaicReview

app = MosaicReview("cubes.lis", cnet_path="control.net")
app.panel()  # renders inline
```

```python
from isistools.apps.tiepoint_review import TiepointReview

app = TiepointReview("cubes.lis", "control.net")
app.panel()
```

### Low-level API

```python
from isistools.io.footprints import load_footprints
from isistools.io.controlnet import load_cnet
from isistools.io.cubes import load_cube
from isistools.plotting.footprint_map import footprint_map
from isistools.plotting.image_viewer import image_plot, image_with_cnet
from isistools.plotting.cnet_overlay import cnet_to_geodataframe

# Load data
gdf = load_footprints("cubes.lis")
cnet = load_cnet("control.net")
da = load_cube("image.cub")

# Plot footprints
footprint_map(gdf)

# Plot image with control points
image_with_cnet(da, cnet, serial_number="MRO/HIRISE/...")

# Convert cnet to GeoDataFrame for map overlay
cnet_gdf = cnet_to_geodataframe(cnet)
```

## Requirements

- Cubes must have `footprintinit` already run for footprint display
- Control networks in ISIS3 binary format (compatible with jigsaw)
- GDAL with ISIS3 driver support (for reading .cub files via rioxarray)

## Dependencies

Core: geopandas, hvplot, holoviews, datashader, panel, rioxarray, plio, pvl, typer

## License

MIT
