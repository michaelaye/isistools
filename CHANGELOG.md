# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.8.0] - 2026-04-12

Body-agnostic refactor. Before 0.8.0, csm2map silently assumed Mars —
a JANUS, LRO, or Europa Clipper cube would have been projected with
Mars radii and the user would have had no way to notice. 0.8.0 pulls
the target body's ellipsoid directly from ALE's ISD and threads it
through the pipeline with zero hardcoded literals.

Output on Mars inputs is **bit-for-bit identical** to 0.7.1. Verified
end-to-end on a full-length F05 CTX cube (557M pixels): `cmp -s` on the
two GeoTIFFs returns success.

### BREAKING

- **`isistools.processing.camera.load_camera()` now returns a tuple**
  `(csmapi.RasterGM, TargetBody)` instead of just the model. Any code
  that unpacks `load_camera(cube)` must be updated. The CLI is
  unaffected. Only `isistools.processing.project.project()` called
  this symbol directly, and it has been updated to match.
- **`isistools.processing.dem.DemRadiusSampler.__init__`'s
  `fallback_radius` is now a required keyword argument** — the old
  Mars-specific default (`3389526.7`) is gone. Callers must pass the
  target body's mean radius explicitly.
- **`isistools.geo.projections.mapping_to_crs()` raises
  `ValueError`** when the Mapping group lacks `EquatorialRadius`.
  Previous versions silently defaulted to Mars radii, which silently
  mis-projected any non-Mars Mapping group. If you have a MAP file
  without explicit radii, add them.
- **`isistools.processing.camera.get_target_radii()` has been
  removed.** Its body-specific information is now carried on the
  `TargetBody` returned from `load_camera()`. There was no public
  caller of `get_target_radii` outside the csm2map pipeline itself.
- **`isistools.processing.grid.grid_from_map_file()` raises
  `ValueError`** when a MAP file uses `Scale` (pixels/degree) but
  lacks `EquatorialRadius`. Scale→resolution conversion needs the
  body's equatorial radius and previously silently used Mars.

### Added

- **`isistools.processing.camera.TargetBody`** (new public
  dataclass): frozen dataclass describing a target body's
  ellipsoid and identity:
    - `name` (str, e.g. `"MARS"`, `"EUROPA"`)
    - `naif_id` (int, e.g. 499, 502, 301)
    - `radius_equatorial_m`, `radius_polar_m`, `radius_mean_m`
      (all in meters)
  Built via the `TargetBody.from_isd(isd_dict, target_name=...)`
  classmethod, which parses ALE's native ISD format and converts
  km → m automatically. Includes a cross-check that validates the
  ISD's top-level `radii` dict against
  `naif_keywords.BODY<code>_RADII` and raises `ValueError` if they
  disagree by more than 1 meter (catches stale SPICE blobs and
  corrupted cubes).
- **csm2map now prints the target body on startup** — the output
  log now shows `Target: <NAME> (NAIF <id>)  radii eq=... polar=...`
  so the user can confirm the right body is being used. No more
  silent Mars assumption.
- **`tests/test_target_body.py`** (12 new regression tests): covers
  Mars / Moon / Europa / a hypothetical `BODY999` through
  `TargetBody.from_isd`, unit conversion (km vs m), the
  ISD-vs-BODY_RADII cross-check, name uppercasing, and the
  `mapping_to_crs` hardening. Test count: 22 → 34.

### Changed

- **`load_camera()` reads and caches the ISD JSON once** — the
  same string is both handed to `csmapi.Isd()` (via the JSON file
  on disk) and parsed to a Python dict (for `TargetBody`). No
  duplicate ALE calls, no second SPICE query.
- **`project.project()` pipeline simplified**: the old separate
  `get_target_radii()` stage is gone. The `TargetBody` comes out
  of `load_camera()` already populated, so the `target_radii`
  timing stage has been removed from `--profile` output.
- **`_build_grid()` constructs its default projection string at
  runtime** from `body.radius_equatorial_m` / `body.radius_polar_m`
  instead of the hardcoded Mars ellipsoid. When no MAP file and no
  explicit `--projection` flag are given, csm2map now picks the
  correct body automatically.

### Fixed

- **Hardcoded Mars radii removed from four sites**:
    - `processing/camera.py::get_target_radii` (deleted entirely)
    - `processing/project.py::_build_grid` default projection string
    - `processing/dem.py::DemRadiusSampler.fallback_radius` default
    - `geo/projections.py::mapping_to_crs` silent default
    - `processing/grid.py::grid_from_map_file` Scale-handling
      branch fallback
  A repo-wide grep for `BODY499`, `3396190`, `3376200`, `3389526`
  in `src/` now returns **zero matches** except for one comment in
  `grid.py` documenting the previous behavior for future readers.

### Validation

- **F05 CTX bit-identical regression**: the 0.8.0 output for the
  F05 full-length CTX cube (1.0 GB input, 557M-pixel output) is
  byte-for-byte identical to the 0.7.1 reference output. Verified
  with `cmp -s`. The body-agnostic refactor has zero behavior
  change on Mars data, which is the only target csm2map has ever
  been run against.

## [0.7.1] - 2026-04-12

Documentation-only patch release. No code behavior changes — all
csm2map pipeline outputs are bit-identical to 0.7.0.

### Documentation

- **Corrected `--clip-to-footprint` documentation.** The flag was
  originally documented as an "ISIS cam2map compatibility mode" that
  matched `cam2map` by clipping csm2map output to the `footprintinit`
  polygon. That framing was based on a working hypothesis later
  disproved empirically: ISIS `cam2map` ignores the polygon entirely
  (stripping the polygon from the cube and re-running `cam2map`
  produces bit-identical output — see `docs/csm2map.qmd § "The
  footprintinit polygon precision story"` for the full narrative).
  All affected documentation has been rewritten to drop the
  cam2map-matching claim and reframe the flag as an escape hatch for
  downstream tooling that explicitly wants a polygon-shaped output
  mask. Affected surfaces:
    - `isistools csm2map --help` — Typer option help text and
      command docstring no longer claim ISIS-compatibility.
    - `docs/csm2map.qmd` — Purpose paragraph, Usage example comment,
      options table row, and the former "Pixel-perfect match"
      paragraph rewritten. Cross-references the empirical-disproof
      section.
    - `docs/csm2map-design.md` §6 Known limitations gains a new
      entry 7 documenting the flag as a historical
      hypothesis-test-reject artifact, flagged as paper-worthy.
    - `src/isistools/processing/project.py` — `clip_to_footprint`
      parameter docstring, the in-function console message, and
      the `_rasterize_footprint` helper docstring all corrected.

### Known limitations (unchanged from 0.7.0)

The `--clip-to-footprint` flag itself is retained in 0.7.1 for
backward compatibility. Deprecation or removal is being considered
for a future release.

## [0.7.0] - 2026-04-12

### Added

- **`csm2map` CLI command**: CSM-based replacement for ISIS `cam2map`.
  Map-projects an ISIS cube into a GeoTIFF using the Community Sensor
  Model (via `ale` + `usgscsm` + `csmapi`) instead of ISIS's CSPICE camera.
  Validated against ISIS 9.0.0 at 99.95% coverage match and 100% agreement
  within 0.01 DN on CTX; runs 5–13× faster than ISIS `cam2map` depending
  on cube length. Optional feature gated behind the `[csm]` extra — base
  isistools users pay no new dependency cost.
- **`csm2map-compare` CLI command**: numeric validation of a csm2map
  GeoTIFF against an ISIS cam2map reference cube.
- **Reads an existing DEM as a shape model during projection**
  (`--shape-model auto|ellipsoid|<path>`): csm2map consumes a DEM cube
  (e.g. the MOLA radius DEM for Mars) to look up per-pixel body radii when
  back-projecting through the CSM sensor. `auto` reads the input cube's
  `Kernels.ShapeModel` and opens the same DEM ISIS `cam2map` would use,
  matching ISIS's default behavior. This is a *read* path — csm2map uses
  a DEM, it does not produce one.
- **Jigsaw-aware SPICE source** (`--spice-source isis|naif|auto`, default
  `isis`): reads SPICE pointing/position from the cube's embedded blobs
  instead of the live NAIF kernels, which is the only correct choice after
  `jigsaw update=true` since jigsaw updates the blobs but NOT the live
  kernels. `naif` and `auto` are available for comparisons against
  pre-jigsaw geometry.
- **Stage-timing profiler** (`--profile`): per-stage wall time breakdown
  (camera load, DEM open, coord transform, read input, resample, write).
- **`docs/csm2map.qmd`**: full user chapter covering purpose, installation
  (including the Apple Silicon `csmapi` build-from-source gotcha),
  usage, options table, pipeline description, and ISIS-compatibility
  validation.
- **`docs/csm2map-design.md`**: design document capturing the resolution
  decisions, validation methodology, and performance trade-offs (for
  citation in an upcoming paper).
- **`scripts/benchmark_csm2map.sh`**: parameterized benchmark harness that
  runs `cam2map` and `csm2map` back-to-back on any CTX cube and emits a
  speedup summary. Uses `pvl.loads()` on `camrange` output instead of
  brittle grep/awk parsing.
- **Processing layer** (`src/isistools/processing/`): new subpackage
  containing the csm2map pipeline — `camera.py` (CSM model loader),
  `grid.py` (output raster grid from MAP file or params), `transform.py`
  (coarse-grid coordinate map with vectorized+threaded bilinear upsample),
  `resample.py` (threaded `scipy.ndimage.map_coordinates`), `writers.py`
  (ZSTD GeoTIFF writer), `project.py` (pipeline orchestrator),
  `dem.py` (lazy windowed DEM radius sampler).
- **PyPI classifiers** in `pyproject.toml`: license, OS, Python versions
  3.10–3.13, Scientific/Astronomy, GIS, Image Processing topics.
- **`[tool.pytest.ini_options]`** with a registered `slow` marker.
- **`tests/test_cli.py`**: smoke tests for every registered CLI command's
  `--help`, plus a regression test for the `overlaps --png` NameError bug
  (item 3 below).

### Changed

- **CLI renamed `cam2map` → `csm2map`** to avoid confusion with the ISIS
  command and to make clear that this is the CSM-based pipeline.
- **Default coarse step 16 → 32** in `compute_transform_coarse`. Validation
  against the dense path shows sub-pixel accuracy at step=32 on CTX,
  roughly halving the CSM-call count with no observable quality cost.
- **ZSTD-compressed tiled GeoTIFF output**, written with multi-threaded
  encoding (`num_threads=ALL_CPUS`).
- **Float32 throughout the resample pipeline** — coordinate maps,
  interpolated coordinates, and output pixels are all `float32`, halving
  the memory bandwidth vs the previous `float64` path.
- **Bilinear coordinate upsample rewritten**: vectorized + thread-striped
  `_bilinear_upsample_pair()` replacing `scipy.ndimage.zoom`, amortizing
  the per-row index math across the line and sample channels. Combined
  with the other CPU optimizations above, this delivered a 3× end-to-end
  speedup vs the 0.6.0 prototype.
- **Threaded resample**: `scipy.ndimage.map_coordinates` releases the GIL,
  so the per-band resample is now split into horizontal stripes processed
  by a `ThreadPoolExecutor`. ~2.5× speedup on the resample stage alone.
- **Fast PVL label parser** (`read_label(fast=True)`, default): reads the
  first 1 MB of the cube and parses it with `pvl.loads()` instead of
  seeking through the whole file. Cut DEM label parsing from 1254 ms to
  15 ms.
- **Removed `knoten` dependency**: the CSM model is now constructed via
  `csmapi.Isd()` + `plugin.constructModelFromISD()` directly, avoiding
  `knoten`'s broken conda packaging.

### Fixed

- **ISIS SIGSEGV under restrictive sandbox**: `campt`/`cam2map`/`camrange`
  and anything else that constructs an ISIS `CubeManager` crashed with
  signal 11 when `RLIMIT_NOFILE` was set to `INT64_MAX`, because
  `CubeManager::p_maxOpenFiles = rlim_cur * 0.60` overflowed to a garbage
  value that corrupted the open-cube cache. `scripts/` and the benchmark
  now set `ulimit -n 4096` up-front. See `scripts/isis_sandbox_fix.md`
  for the full diagnosis.
- **`overlaps --png` crashed with `NameError`** (`cli.py:395`):
  `geopandas` (aliased as `gpd`) was referenced inside the PNG plot
  branch but never imported. Ships in 0.6.0 and earlier. Fixed by adding
  `import geopandas as gpd` next to the lazy matplotlib imports.
  Regression test in `tests/test_cli.py` exercises the PNG branch with
  mocked ISIS/geopandas data and asserts exit code 0.
- **Empty `archive` and `clock_lookup` dead variables** removed from
  `io/cubes.py:get_serial_number()` and
  `apps/mosaic_review.py:_on_image_selected()`. Both were leftovers from
  superseded code paths; the enclosing functions already worked
  correctly without them.
- **Ambiguous variable name `l`** in `plotting/cnet_overlay.py` renamed
  to `line` for readability.

### Documentation

- **Research note** (`scripts/csm_research.md`): CSM governance, the
  verification methodology behind the NGA standard, and a roadmap for
  adding JANUS (ESA JUICE) to the supported sensor set.
- **Polygon-dependence investigation** (`docs/disagreement_analysis.md`):
  empirically verified that ISIS `cam2map` does NOT use the footprint
  polygon stored by `footprintinit`. Two experiments (stripped polygon
  and tight polygon) produced bit-identical outputs, refuting the
  initial hypothesis that the ~870K-pixel coverage gap was polygon-based.
- **Residual disagreement analysis**: the remaining ~18K-pixel difference
  between csm2map and ISIS `cam2map` (after DEM integration) is traced to
  a CSM-vs-CSPICE camera-model floor, not a csm2map bug. Documented in
  `docs/csm2map-design.md` §6.
- **CLAUDE.md**: documents the `py312` / `isis` two-env setup, the
  processing layer architecture, and the csm2map command.

### Known limitations

- `get_target_radii()` currently hardcodes Mars NAIF ID 499. Needs
  generalization for non-Mars targets before JANUS / non-Mars support.
- `csm2map-compare` can produce a 1-row shape mismatch on very long CTX
  strips (e.g. F05, 52,212 vs 52,207 rows) from differing half-pixel
  handling in the lat/lon → pixel rounding. Affects only the comparison
  tool, not the projected output itself.

## [0.6.0] - 2026-03-24

### Added

- **`spiceinit` CLI command**: batch-run ISIS `spiceinit` on all cubes in a
  list file with parallel execution (`-j` flag, default 4 workers). Web kernel
  retrieval enabled by default (`--web`); disable with `-W` / `--no-web`.
- **`overlaps` CLI command**: runs ISIS `findimageoverlaps` and parses the
  output WKB polygons into a GeoDataFrame. Prints a summary table of all
  overlap zones with types and areas. Supports `--png` for quick visualization
  and `--gpkg` for GeoPackage export (loadable in QGIS or notebooks).
- New `isistools.io.overlaps` module with `parse_overlap_list()` function for
  programmatic access to findimageoverlaps output as a GeoDataFrame.

## [0.5.3] - 2026-03-11

### Fixed

- HoloViews shared-axis linking between map and image plots caused axis
  conflicts — added `linked_axes=False`, `shared_axes=False`, and renamed
  image dims from x/y to sample/line.

### Changed

- Cnet point styles: smaller markers (size 3), cross markers for
  registered/unregistered, red circle for ignored points (QA flag).

## [0.5.2] - 2026-03-10

### Changed

- Streamlined README: `pip install isistools` as primary install, detailed
  API docs moved to documentation site, added PyPI badge and docs link.

## [0.5.1] - 2026-03-10

### Fixed

- Inline code invisible in dark mode (superhero theme) — added CSS override.
- README rendered Quarto YAML frontmatter as raw text on GitHub — removed it.

### Added

- Zenodo DOI badge in README.
- Sandstone (light) / Superhero (dark) theme switching in docs.

## [0.5.0] - 2026-03-10

### Added

- **Quarto book documentation** with one page per CLI subcommand, Python API
  reference, and real example figures from CTX test data.
- **GitHub Actions workflow** for building docs and deploying to GitHub Pages.
- **`[project.urls]`** in pyproject.toml for PyPI sidebar links (Homepage,
  Repository, Changelog, Issues).

## [0.4.0] - 2026-03-10

### Added

- **`footprintinit` CLI command**: batch-run ISIS `footprintinit` on all cubes
  in a list file with parallel execution (`-j` flag, default 4 workers).
- **Static PNG export for footprints** (`--png` flag): publication-ready
  footprint overview images via headless matplotlib (Agg backend). Defaults to
  `footprints_overview.png`; override with `--png-path`. Supports `--dpi` and
  `--title` options.
- New `footprint_png()` function in `plotting.footprint_mpl` for programmatic
  PNG export.

### Changed

- Refactored `footprint_mpl.py`: extracted shared `_plot_footprints()` core
  used by both the interactive window viewer and the new PNG exporter.
  `mplcursors` and `QtAgg` backend are now only imported in the interactive
  `footprint_window()` path.
- `footprints` CLI command accepts `--title` option (defaults to cubelist stem).

## [0.3.0] - 2025-02-12

### Added

- **Disk caching via `diskcache`**: repeated runs skip expensive I/O when
  source files are unchanged. Cached items (at `~/.cache/isistools/`):
  - Per-cube footprint records (PVL parsing + polygon blob reading)
  - Control network DataFrames (protobuf decoding + status classification)
  - Per-cube campt coordinate conversions (subprocess calls)
- Cache keys include file `mtime_ns` so entries auto-invalidate on
  source file changes.
- New module `isistools.io.cache` with `get_cache()` accessor.
- New dependency: `diskcache>=5.6`.

### Changed

- `load_cnet()` now returns a plain `pd.DataFrame` (stripped of plio's
  `IsisControlNetwork` protobuf internals) for cache compatibility and
  reduced memory footprint.

## [0.2.2] - 2025-02-12

### Changed

- **Eliminated triple PVL parsing**: cube labels are now parsed once in
  `load_footprints()` and reused for footprint extraction and serial lookup.
  `read_footprint()` accepts an optional pre-parsed `label` parameter.
  `load_footprints()` now extracts `SpacecraftClockCount` during the single
  parse, stored in the `clock` column of the GeoDataFrame.
- **Parallelized campt calls**: `_lonlat_from_campt()` now runs all per-cube
  campt subprocess calls concurrently via `ThreadPoolExecutor`. With N cubes,
  coordinate conversion is roughly N times faster.
- `cnet_to_geodataframe()` accepts optional `clock_lookup` dict to skip
  rebuilding the serial-number-to-cube mapping from labels.

## [0.2.1] - 2025-02-12

### Changed

- Browser footprint map simplified to single hvplot call with `c="short_pid"`
  color mapping instead of per-filename overlay loop.
- Added `ctx_short_pid()` utility in `styles.py`, used by both browser and
  `--win` code paths for consistent 18-char CTX product ID display.
- Registered control points now use black crosses in browser map (matching
  `--win` path); hover disabled on cnet points.
- Plot enlarged to 1200x800 with larger axis labels (20pt), tick labels (14pt),
  and legend text (15pt). Legend below plot in 3 columns with spacing.
- Removed `geo=True` from footprint map; using `data_aspect=1` with free-form
  `BoxZoomTool` for undistorted footprints with flexible zoom.
- Cleaned up footprint hover: removed instrument, target, and filename fields;
  shortened start_time to seconds.

### Fixed

- `FileNotFoundError` on `import panel` in non-`--win` path: added `_ensure_cwd()`
  guard in CLI before any panel/holoviews import.

## [0.2.0] - 2025-02-12

### Added

- **Native matplotlib footprint viewer** (`--win` flag): lightweight alternative
  to the browser-based Panel/Bokeh viewer. Launches a native window with
  zoom/pan toolbar, hover tooltips showing CTX product ID, and a right-side
  legend with per-filename colors.
- **Precise campt coordinate conversion**: pre-jigsaw control networks without
  ground coordinates are converted from sample/line to lon/lat using ISIS
  `campt` in batch mode via kalasiris, providing camera-model-accurate positions.
- Control point overlay on matplotlib footprint viewer: registered points shown
  as black crosses, unregistered as red circles, ignored as gray circles.
- Window centering support for Qt and Tk matplotlib backends.
- New dependencies: `matplotlib`, `mplcursors`, `kalasiris`, `PyQt6`.

### Changed

- Serial lookup utilities (`build_serial_lookup`, `match_serials_to_cubes`)
  moved from `apps.tiepoint_review` to `io.cubes` for reuse without heavy
  imports.
- Holoviews/hvplot/panel imports made lazy in `cnet_overlay.py` — the `--win`
  path no longer loads browser-plotting machinery, fixing startup crashes and
  improving load time.
- `cnet_to_geodataframe()` now accepts optional `cube_paths` parameter for
  campt-based coordinate conversion on level-1 networks.

### Fixed

- `NameError` in `mosaic_review.py` where `_build_serial_lookup` was called but
  never imported.
- `FileNotFoundError` on startup caused by `os.getcwd()` in `param` triggered
  via transitive holoviews imports in the matplotlib code path.

## [0.1.1] - 2025-02-11

### Fixed

- Control network loading: renamed plio columns (`id`→`pointId`,
  `sampleResidual`→`residualSample`, `lineResidual`→`residualLine`) to
  match isistools conventions.
- Serial number matching: replaced broken filename-based heuristic with
  `SpacecraftClockCount` lookup from cube labels, fixing "Could not find
  cube files for this pair" in tiepoint review.
- Duplicate toolbar buttons (wheel_zoom, pan, reset) in image viewer
  resolved via a Bokeh post-render deduplication hook.

## [0.1.0] - 2025-02-11

### Added

- Initial release.
- I/O layer: cube loading (`load_cube`), footprint reading
  (`read_footprint`, `load_footprints`), control network I/O
  (`load_cnet`, `save_cnet`) via plio.
- Plotting layer: rasterized image viewer with percentile contrast
  stretch, footprint map with per-filename Category20 colors and
  clickable mute legend, control network point overlays in both image
  and map space.
- App layer: `MosaicReview` (Qmos replacement) and `TiepointReview`
  (Qnet replacement) as Panel apps.
- Typer CLI with commands: `mosaic`, `tiepoints`, `footprints`,
  `cnet-info`.

[Unreleased]: https://github.com/michaelaye/isistools/compare/v0.8.0...HEAD
[0.8.0]: https://github.com/michaelaye/isistools/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/michaelaye/isistools/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/michaelaye/isistools/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/michaelaye/isistools/compare/v0.5.3...v0.6.0
[0.5.3]: https://github.com/michaelaye/isistools/compare/v0.5.2...v0.5.3
[0.5.2]: https://github.com/michaelaye/isistools/compare/v0.5.1...v0.5.2
[0.5.1]: https://github.com/michaelaye/isistools/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/michaelaye/isistools/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/michaelaye/isistools/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/michaelaye/isistools/compare/v0.2.2...v0.3.0
[0.2.2]: https://github.com/michaelaye/isistools/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/michaelaye/isistools/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/michaelaye/isistools/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/michaelaye/isistools/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/michaelaye/isistools/releases/tag/v0.1.0
