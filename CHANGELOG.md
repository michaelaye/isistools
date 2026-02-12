# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/michaelaye/isistools/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/michaelaye/isistools/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/michaelaye/isistools/releases/tag/v0.1.0
