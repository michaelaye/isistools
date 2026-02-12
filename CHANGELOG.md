# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/michaelaye/isistools/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/michaelaye/isistools/compare/v0.2.2...v0.3.0
[0.2.2]: https://github.com/michaelaye/isistools/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/michaelaye/isistools/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/michaelaye/isistools/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/michaelaye/isistools/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/michaelaye/isistools/releases/tag/v0.1.0
