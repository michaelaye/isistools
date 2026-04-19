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

# csm2map - standalone CLI, CSM-based replacement for ISIS cam2map
# Runs from the py312 conda env (not the isis env): csmapi is built from
# source for osx-arm64 and installed there alongside ale/usgscsm/isistools.
# The isis env is used only for the ISIS reference tools (cam2map, camrange,
# catlab) that our benchmarks/comparisons shell out to.
csm2map input.cub output.tif                    # auto everything
csm2map input.cub output.tif --map equi.map     # use an ISIS MAP file
csm2map input.cub output.tif -r 6.0             # explicit resolution
csm2map compare isis_output.cub csm_output.tif  # validate vs ISIS

# ctxpipe - standalone CLI, Python replacement for ISIS mroctx2isis+ctxcal+ctxevenodd
ctxpipe B04_011267_0983_XN_81S063W.IMG calibrated.tif  # full pipeline
ctxpipe B04_011267_0983_XN_81S063W.IMG out.tif --iof --sun-distance 2.28e8  # I/F units

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

### csm2map subpackage (`csm2map/`)
Self-contained CSM-based replacement for ISIS `cam2map`, with its own standalone CLI (`csm2map` at the shell) and Python API (`from isistools.csm2map import csm2map`). Could be extracted as a standalone package in the future — shared deps are just `io/cubes.py` (`read_label`, `read_isis_cube_raw`) and `io/footprints.py` (`read_footprint`).

Pipeline: camera model → output grid → coordinate transform → resample → GeoTIFF + PVL sidecar.

- `camera.py` — `load_camera()` returns `(csmapi.RasterGM, TargetBody)`. `TargetBody` frozen dataclass carries body ellipsoid + NAIF ID from ALE ISD. `compute_ground_sample_distance()` for auto-resolution. `ground_to_image_batch()` for coordinate mapping.
- `grid.py` — `OutputGrid` dataclass defining the output raster (CRS, affine, dimensions). `grid_from_map_file()` handles all ISIS conventions (planetographic, positive-west, 180/360 domain). `grid_from_params()` with ISIS-compatible snap rule.
- `transform.py` — `CoordinateMap` mapping output → input pixels. Coarse-grid + vectorized bilinear interpolation (production, step=32) or dense (validation).
- `resample.py` — `scipy.ndimage.map_coordinates` with threaded stripe processing. Supports nearest/bilinear/bicubic.
- `writers.py` — GeoTIFF (ZSTD, tiled) + ISIS-compatible Mapping PVL sidecar.
- `pipeline.py` — `csm2map()` function orchestrating the full pipeline. Auto-resolution from camera GSD, auto-bounds via circular-statistics footprint derivation, auto-centered projection.
- `projections.py` — `mapping_to_crs()` converts ISIS Mapping group to `pyproj.CRS`. `planetographic_to_planetocentric()` and longitude normalization helpers.
- `compare.py` — Pixel-level comparison of csm2map output vs ISIS cam2map reference.
- `cli.py` — Standalone Typer app with `project` and `compare` commands.

### ctxpipe subpackage (`ctxpipe/`)
Pure-Python replacement for the ISIS CTX calibration chain (mroctx2isis + ctxcal + ctxevenodd). Reads PDS3 EDRs directly, no ISIS installation needed. Validated pixel-exact against ISIS output (max relative error < 1e-6).

- `ingest.py` — `ingest_ctx_edr()` reads PDS3 EDR, applies SQROOT decompression (8-bit → 12-bit LUT), extracts dark pixels. Returns `(image, CTXMetadata)`.
- `calibrate.py` — `calibrate()` applies dark current subtraction (per-channel A/B for summing=1), flat-field correction from `$ISISDATA/mro/calibration/ctxFlat_NNNN.cub`, optional I/F conversion. Constants: w0=3660.5, perihelion=2.07e8 km.
- `evenodd.py` — `correct_evenodd()` fixes alternating-column striping. Only for summing=1. Correction = half the even-odd mean difference, added to odd / subtracted from even.
- `pipeline.py` — `ctxpipe()` orchestrates ingest → calibrate → evenodd. Optional GeoTIFF output.
- `cli.py` — Standalone Typer app with `calibrate` command. Entry point: `ctxpipe`.

### hirisepipe subpackage (`hirisepipe/`)
Pure-Python replacement for the ISIS HiRISE RED calibration chain (hical + histitch + cubenorm). Reads ISIS cubes (from hi2isis), no ISIS needed for calibration. Validated against ISIS hical to 0.08% mean relative error.

- `hical.py` — `hical()` applies the full 10-module calibration chain: ZeroBufferSmooth, ZeroBufferFit, ZeroReverse, ZeroDark, GainLineDrift, GainNonLinearity, GainChannelNormalize, GainFlatField, GainTemperature, GainUnitConversion. Two-pass algorithm: first pass subtracts drift/dark/reverse, second pass applies gain/flat/temp using line median for non-linearity. Reads calibration matrices from `$ISISDATA/mro/calibration/matrices/`.
- `stitch.py` — `stitch_channels()` stitches CCD channels 0 and 1 with optional balance correction (seam-edge average ratio).
- `cubenorm.py` — `cubenorm()` normalizes column-to-column variations by dividing each column by its median.
- `pipeline.py` — `process_ccd()` orchestrates hical → histitch → cubenorm for a single CCD.

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
