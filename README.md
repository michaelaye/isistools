# isistools

[![DOI](https://zenodo.org/badge/1155816216.svg)](https://doi.org/10.5281/zenodo.18944803)
[![PyPI](https://img.shields.io/pypi/v/isistools)](https://pypi.org/project/isistools/)

Python-based review tools for ISIS3 coregistration workflows. Replaces Qmos
and Qnet with modern, interactive visualization using HoloViews, Panel, and
datashader.

**[Documentation](https://michaelaye.github.io/isistools/)**

## Installation

```bash
pip install isistools
```

For development:

```bash
git clone https://github.com/michaelaye/isistools.git
cd isistools
pip install -e ".[dev]"
```

## Quick start

```bash
# Interactive footprint map in the browser
isistools footprints cubes.lis

# Export publication-ready PNG
isistools footprints cubes.lis --png --title "My Mosaic"

# Mosaic review with control network (Qmos replacement)
isistools mosaic cubes.lis --cnet control.net

# Tie point review (Qnet replacement)
isistools tiepoints cubes.lis control.net

# Batch footprintinit (parallel)
isistools footprintinit cubes.lis -j 8

# Control network summary
isistools cnet-info control.net
```

## Why?

ISIS's Qmos and Qnet have several pain points:

- **Qmos** requires level-1 images but displays footprints in map projection,
  and renders image content very slowly.
- **Qnet** requires level-2 images but displays them flipped in detector
  readout order.
- **Both** use a color scheme where registered control points are nearly
  invisible while unregistered points are prominent.

isistools fixes all of these. See the
[documentation](https://michaelaye.github.io/isistools/) for details on each
command, the Python API, and example figures.

## License

MIT
