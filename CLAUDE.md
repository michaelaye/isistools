# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

isistools replaces ISIS3's Qmos and Qnet with Python-based interactive review tools for planetary image coregistration workflows. It reads ISIS .cub files and binary control networks, visualizing footprints, images, and tie points using HoloViews/Panel/Bokeh served in a browser. It also provides `csm2map`, a CSM-based map-projection tool that replaces ISIS `cam2map`.

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

# csm2map - CSM-based replacement for ISIS cam2map
# Runs from the py312 conda env (not the isis env): csmapi is built from
# source for osx-arm64 and installed there alongside ale/usgscsm/isistools.
# The isis env is used only for the ISIS reference tools (cam2map, camrange,
# catlab) that our benchmarks/comparisons shell out to.
#   /Users/maye/miniforge3/envs/py312/bin/isistools csm2map ...
#   conda activate isis && cam2map from=... to=...
isistools csm2map input.cub output.tif --map equi.map
isistools csm2map input.cub output.tif -r 6.0
isistools csm2map input.cub output.tif --map equi.map --dense --validate
isistools csm2map-compare isis_output.cub csm_output.tif

# Quarto docs (front matter in README.md, no _quarto.yml needed)
quarto render README.md
quarto preview README.md
```

## Architecture

The codebase follows a three-layer pattern: **I/O → Plotting → Apps**

### I/O layer (`io/`)
- `footprints.py` — Reads footprint polygons from ISIS cubes. The polygon is stored as a WKT/GML blob at a byte offset specified by the `Object = Polygon` block in the PVL label (`StartByte`/`Bytes` fields). **Not** a `^Polygon` pointer.
- `controlnet.py` — Wraps plio's `from_isis()`/`to_isis()` (not `read_network`/`write_network`, those don't exist). Adds `residual_magnitude` and `status` columns.
- `cubes.py` — Loads .cub files as xarray DataArrays via rioxarray/GDAL. Normalizes orientation to north-up for level-2 cubes. Also provides `read_isis_cube_raw()` for direct binary reading (used by csm2map pipeline) and ISIS special pixel constants.

### Plotting layer (`plotting/`)
- `footprint_map.py` — Per-filename colored polygons with Category20 palette, clickable mute legend via `legend_opts`. Each filename is a separate hvplot overlay.
- `image_viewer.py` — Datashader-rasterized image display with percentile contrast stretch.
- `cnet_overlay.py` — Control points in both image space (sample/line) and map space (lon/lat), with status-based coloring (registered=green, unregistered=red, ignored=gray).
- `styles.py` — All visual constants: point colors/sizes by status, footprint colors.

### App layer (`apps/`)
- `mosaic_review.py` — MosaicReview class: footprint map + image browser + cnet overlay. Qmos replacement.
- `tiepoint_review.py` — TiepointReview class: side-by-side image pairs with shared tie points. Qnet replacement.
- `components.py` — Reusable Panel widgets (cube list selector, cnet selector, info panels).

### Processing layer (`processing/`)
CSM-based replacement for ISIS `cam2map`, exposed as the `csm2map` CLI command. Pipeline: camera model → output grid → coordinate transform → resample → GeoTIFF.
- `camera.py` — Loads CSM RasterGM sensor models from spiceinit'd cubes via ale/usgscsm (no knoten dependency). `load_camera()` returns `(csmapi.RasterGM, TargetBody)` — the `TargetBody` frozen dataclass carries the target body's ellipsoid (radii, NAIF ID, mean radius) parsed from ALE's ISD. Provides `ground_to_image_batch()` for coordinate mapping (the performance bottleneck — Python loop over CSM calls).
- `grid.py` — `OutputGrid` dataclass defining the output raster (CRS, affine transform, dimensions). Built from ISIS MAP files or explicit parameters.
- `transform.py` — `CoordinateMap` mapping output pixels to input pixels. Two strategies: dense (every pixel, slow) and coarse-grid + bilinear interpolation (production path, step=16 default).
- `resample.py` — Applies coordinate map via `scipy.ndimage.map_coordinates`. Supports nearest/bilinear/bicubic.
- `writers.py` — GeoTIFF output via rasterio.
- `project.py` — Orchestrates the full pipeline. Also contains `_derive_ground_range()` for auto-computing bounds from camera model.

### Geo layer (`geo/`)
- `projections.py` — `mapping_to_crs()` converts ISIS Mapping group to `pyproj.CRS`. `mapping_to_proj4()` is a thin wrapper. `ISIS_TO_PROJ4` maps 12 ISIS projection names to proj4 IDs.

### CLI (`cli.py`)
Typer-based. Each command constructs an app object and calls `.serve()`. Entry point: `isistools = "isistools.cli:app"`. The `csm2map` command guards its CSM imports with try/except for a clear error if `[csm]` extras are not installed.

## Versioning and Releases

Single source of truth: `__init__.py` defines `__version__`. `pyproject.toml` reads it dynamically via `[tool.hatch.version]`. **Only edit `__init__.py` when bumping versions.**

When releasing a new version, always do all three:
1. Bump `__version__` in `src/isistools/__init__.py`
2. Update `CHANGELOG.md` following [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) format (move Unreleased items to the new version section, add comparison links)
3. Tag the commit as `v{version}`

## Key Dependencies and Gotchas

- **plio**: API is `from_isis(path)` and `to_isis(df, path)`. Older names like `read_network`/`write_network` do not exist. Plio columns are renamed in `load_cnet()`: `id`→`pointId`, `sampleResidual`→`residualSample`, `lineResidual`→`residualLine`.
- **Serial numbers**: ISIS serial numbers (e.g., `MRO/CTX/0910464726:234`) use spacecraft clock counts, not product IDs. Match cubes via `SpacecraftClockCount` from the label's Instrument group.
- **pvl**: ISIS labels are parsed as `PVLModule`. The `Polygon` object has `StartByte` (1-based) and `Bytes`.
- **hvplot**: `active_tools` and `default_tools` are not valid parameters for hvplot plots — they produce warnings. Use Bokeh hooks or `legend_opts` instead.
- **hvplot `by=` with `cmap=`**: Does not reliably apply fill/line colors to polygons. Use explicit per-group overlays with `fill_color`/`line_color` instead.
- **Bokeh legend**: Use `legend_opts={"click_policy": "mute"}` on `hv.opts.Overlay` for clickable legends. No Bokeh hooks needed.
- **ISIS cubes**: Require GDAL with ISIS3 driver support. Cubes must have `footprintinit` run before footprints can be read.
- **rioxarray**: Import `rioxarray` to register the `.rio` accessor, even if not called directly.
- **ale**: Use `ale.loads(cube_path)` (not `ale.load()` which is deprecated). Returns ISD JSON string. Install from conda-forge.
- **usgscsm/csmapi**: We call the CSM plugin API directly (no knoten dependency — it has broken conda packaging). `import usgscsm` registers the plugin, then `csmapi.Isd()` + `plugin.constructModelFromISD()` builds the model. Install via `conda install -c conda-forge usgscsm`.
- **csmapi**: Key classes: `RasterGM`, `ImageCoord(line, samp)`, `EcefCoord(x, y, z)`. CSM is 0-based (pixel center at 0.0, 0.0); ISIS is 1-based (pixel center at 1, 1). **Half-pixel offset matters.**
- **Coordinate conventions (csm2map)**: CSM and pyproj use planetocentric lat, positive-east lon. ISIS MAP files may specify planetographic lat or positive-west lon — 0.8.1 adds conversion via `planetographic_to_planetocentric()` and `normalize_longitude()` in `geo/projections.py`. `grid_from_map_file` reads `LatitudeType`, `LongitudeDirection`, and `LongitudeDomain` from the MAP file and normalizes automatically, emitting a `UserWarning` when conversion happens.
- **Target body handling**: 0.8.0+ is body-agnostic. `TargetBody.from_isd()` extracts radii and NAIF ID from ALE's ISD; no Mars-specific constants remain in `src/`. Before 0.8.0, `get_target_radii()` hardcoded Mars NAIF ID 499.

## Style

- Line length: 99 (configured in `[tool.ruff]`)
- Lint rules: E, F, I, W
- Type hints on all public functions
- src layout: package lives under `src/isistools/`
