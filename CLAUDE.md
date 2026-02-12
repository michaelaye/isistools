# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

isistools replaces ISIS3's Qmos and Qnet with Python-based interactive review tools for planetary image coregistration workflows. It reads ISIS .cub files and binary control networks, visualizing footprints, images, and tie points using HoloViews/Panel/Bokeh served in a browser.

## Commands

```bash
pip install -e .                    # Install in development mode
pytest tests/                       # Run all tests
pytest tests/test_io.py::TestClassifyPointStatus::test_ignored_point  # Single test
ruff check src/                     # Lint
ruff format src/                    # Format

# CLI (requires ISIS cubes with footprintinit already run)
isistools mosaic cubes.lis --cnet control.net
isistools tiepoints cubes.lis control.net
isistools footprints cubes.lis
isistools cnet-info control.net

# Quarto docs (front matter in README.md, no _quarto.yml needed)
quarto render README.md
quarto preview README.md
```

## Architecture

The codebase follows a three-layer pattern: **I/O → Plotting → Apps**

### I/O layer (`io/`)
- `footprints.py` — Reads footprint polygons from ISIS cubes. The polygon is stored as a WKT/GML blob at a byte offset specified by the `Object = Polygon` block in the PVL label (`StartByte`/`Bytes` fields). **Not** a `^Polygon` pointer.
- `controlnet.py` — Wraps plio's `from_isis()`/`to_isis()` (not `read_network`/`write_network`, those don't exist). Adds `residual_magnitude` and `status` columns.
- `cubes.py` — Loads .cub files as xarray DataArrays via rioxarray/GDAL. Normalizes orientation to north-up for level-2 cubes.

### Plotting layer (`plotting/`)
- `footprint_map.py` — Per-filename colored polygons with Category20 palette, clickable mute legend via `legend_opts`. Each filename is a separate hvplot overlay.
- `image_viewer.py` — Datashader-rasterized image display with percentile contrast stretch.
- `cnet_overlay.py` — Control points in both image space (sample/line) and map space (lon/lat), with status-based coloring (registered=green, unregistered=red, ignored=gray).
- `styles.py` — All visual constants: point colors/sizes by status, footprint colors.

### App layer (`apps/`)
- `mosaic_review.py` — MosaicReview class: footprint map + image browser + cnet overlay. Qmos replacement.
- `tiepoint_review.py` — TiepointReview class: side-by-side image pairs with shared tie points. Qnet replacement.
- `components.py` — Reusable Panel widgets (cube list selector, cnet selector, info panels).

### CLI (`cli.py`)
Typer-based. Each command constructs an app object and calls `.serve()`. Entry point: `isistools = "isistools.cli:app"`.

## Versioning

Single source of truth: `__init__.py` defines `__version__`. `pyproject.toml` reads it dynamically via `[tool.hatch.version]`. **Only edit `__init__.py` when bumping versions.**

## Key Dependencies and Gotchas

- **plio**: API is `from_isis(path)` and `to_isis(df, path)`. Older names like `read_network`/`write_network` do not exist. Plio columns are renamed in `load_cnet()`: `id`→`pointId`, `sampleResidual`→`residualSample`, `lineResidual`→`residualLine`.
- **Serial numbers**: ISIS serial numbers (e.g., `MRO/CTX/0910464726:234`) use spacecraft clock counts, not product IDs. Match cubes via `SpacecraftClockCount` from the label's Instrument group.
- **pvl**: ISIS labels are parsed as `PVLModule`. The `Polygon` object has `StartByte` (1-based) and `Bytes`.
- **hvplot**: `active_tools` and `default_tools` are not valid parameters for hvplot plots — they produce warnings. Use Bokeh hooks or `legend_opts` instead.
- **hvplot `by=` with `cmap=`**: Does not reliably apply fill/line colors to polygons. Use explicit per-group overlays with `fill_color`/`line_color` instead.
- **Bokeh legend**: Use `legend_opts={"click_policy": "mute"}` on `hv.opts.Overlay` for clickable legends. No Bokeh hooks needed.
- **ISIS cubes**: Require GDAL with ISIS3 driver support. Cubes must have `footprintinit` run before footprints can be read.
- **rioxarray**: Import `rioxarray` to register the `.rio` accessor, even if not called directly.

## Style

- Line length: 99 (configured in `[tool.ruff]`)
- Lint rules: E, F, I, W
- Uses `from __future__ import annotations` throughout
- Type hints on all public functions
- src layout: package lives under `src/isistools/`
