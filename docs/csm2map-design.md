# csm2map: design document

This document captures the design decisions behind `csm2map`, the
CSM-based ISIS-`cam2map` replacement in `isistools`. Its intended
audience is (a) a reviewer of the companion paper, (b) a future
maintainer, and (c) a user who wants to understand *why* the tool
behaves the way it does beyond the user-facing chapter in the Quarto
book.

Document status: live. Updated alongside the code. If a decision
changes, this document changes with it.

## 1. Purpose and scope

`csm2map` projects a spiceinit'd ISIS cube into a map-projected
GeoTIFF using the Community Sensor Model (CSM) API via `ale` and
`usgscsm`, instead of ISIS's CSPICE-based `Isis::Camera` class. It is
functionally equivalent to `cam2map` for the core map-projection case
and validated at pixel level against ISIS 9.0.0 output on multiple
MRO CTX test cases.

The tool exists for two reasons:

1. **Ergonomics and speed.** A pure Python + numpy + scipy pipeline
   can be faster than ISIS `cam2map` (~23× on our MRO CTX benchmark)
   because it can skip ISIS's `ProcessRubberSheet` quad-tree and use
   modern vectorized scientific Python. It also composes naturally
   with other Python tools: rasterio, rioxarray, GDAL, dask,
   scikit-image, etc.
2. **Portability to missions ISIS doesn't support.** Writing a new
   ISIS camera class (C++, C++ review cycle, ISIS release cycle) is
   expensive. ALE drivers (Python, ~150-300 LOC per instrument) are
   much cheaper. A CSM-native projection tool is the natural front end
   for any mission where the CSM + ALE path is already the fastest
   route to a working camera model.

The scope is **map projection only**: input ISIS cube → output
map-projected GeoTIFF. It is explicitly not a mosaicking tool, not a
bundle adjuster, not a radiometric calibrator, and not a stereo
engine. Those are separate tools that can consume `csm2map` output
or feed it different cameras.

## 2. Architecture

The pipeline has seven stages, each measurable independently via the
`--profile` flag:

```
  ┌──────────────────┐
  │  input ISIS cube │
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ load_camera      │   ALE → ISD JSON → usgscsm plugin → RasterGM
  │ (0.6 s)          │
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ target_radii +   │   Reads body radii from cube's NaifKeywords
  │ dem_open         │   and lazily opens the shape-model DEM
  │ (0.2 s)          │
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ build_grid       │   Parses MAP PVL, projects lat/lon corners,
  │ (0.01 s)         │   computes pixel-aligned affine
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ coord_transform  │   For each coarse-grid output pixel:
  │ (0.6 s)          │     inverse-project to lat/lon
  │                  │     look up DEM radius
  │                  │     call CSM groundToImage
  │                  │   Then threaded-bilinear upsample to full res
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ read_input       │   rasterio.read() or np.fromfile (Tile vs BSQ)
  │ (0.15 s)         │
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ resample         │   scipy.ndimage.map_coordinates, threaded
  │ (1.0 s)          │   across horizontal stripes of the output
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ write_output     │   rasterio GeoTIFF, ZSTD compression,
  │ (0.13 s)         │   tiled 256×256, multi-threaded encoder
  └────────┬─────────┘
           │
           ▼
  ┌──────────────────┐
  │ map-projected    │
  │ GeoTIFF          │
  └──────────────────┘
```

Total: ~2.8 s per MRO CTX image on a recent Apple Silicon laptop,
compared to ~65 s for ISIS `cam2map` on the same input.

All seven stages are pure-Python-orchestrated; the heavy lifting is in
ALE (Python reading SPICE via CSPICE), `csmapi.RasterGM` (C++ CSM),
`scipy.ndimage.map_coordinates` (C), `pyproj` (C PROJ), and `rasterio`
(C GDAL).

## 3. Design decisions

### 3.1 CSM camera model instead of ISIS `Isis::Camera`

Decision: use `usgscsm::UsgsAstroLineScanSensorModel` (or framing
model, depending on the sensor) rather than ISIS's CSPICE-based
`Isis::Camera` subclasses.

Rationale:

- The CSM path is fully Python-accessible (`csmapi` via SWIG) without
  linking `libisis`.
- CSM is the emerging interoperability standard across planetary
  photogrammetry tools (USGS, ASP, SOCET GXP).
- Published validation (Laura et al., *Earth and Space Science* 2020,
  doi:10.1029/2019EA000713; Laura et al., *Remote Sensing* 16(4):648,
  2024) confirms CSM and ISIS Camera agree to sub-pixel accuracy on
  the instruments we care about (CTX, HiRISE, HRSC, LROC NAC).
- Residual disagreement at the image boundaries (~14-18 K pixels out
  of 46 M on our CTX test cube) is a structural difference between
  `UsgsAstroLineScanSensorModel` and ISIS `CTXCamera` at the
  line-scan time-domain edges, not a bug in either. It manifests as
  sub-pixel differences in where each implementation converges when
  solving "at what time was the spacecraft pointing at this ground
  point?". Documented and accepted as a noise floor.

Consequence: csm2map DN output agrees with ISIS `cam2map` output to
within 0.01 DN on 100% of overlapping pixels, std 0.0011, across two
different CTX cubes. The remaining ~0.04% coverage difference is the
camera-model boundary noise floor.

### 3.2 ALE with `only_isis_spice=True` for jigsaw awareness

Decision: we pass `only_isis_spice=True` to `ale.loads()` by default,
forcing ALE to read SPICE pointing from the cube's embedded `Table`
blobs (the `IsisSpice` driver path) rather than from live NAIF
kernels.

Rationale: ALE's default driver-selection heuristic (`sort_drivers` in
`ale/drivers/__init__.py`) prefers `NaifSpice` drivers — *"drivers
that use external ephemeris data are tested before drivers that use
attached ephemeris data"*. For any pipeline that runs `jigsaw
update=true`, this default is actively wrong: jigsaw updates the
cube's embedded blobs but does NOT rewrite the live NAIF kernels.
ALE's default behavior would silently discard jigsaw's bundle
adjustment results and project with pre-jigsaw geometry. ALE itself
provides `only_isis_spice=True` for exactly this case; we just use it
by default.

Exposed as `--spice-source {isis,naif,auto}` with `isis` as the
default. The flag is documented in both the Quarto chapter and the
CLI `--help`.

Consequence: csm2map produces correct output for jigsaw-adjusted
cubes out of the box, without requiring users to know about this ALE
heuristic gotcha.

### 3.3 DEM shape model, auto-resolved from the cube label

Decision: by default (`--shape-model auto`), we read the cube's
`Kernels.ShapeModel` keyword, resolve `$ISISDATA` and `$base` path
variables against the environment, and use the referenced DEM cube as
the surface-radius source via a windowed rasterio read. If no
shape model is specified (or if the DEM file is missing), we fall
back to a spherical surface with radius equal to the mean of the
cube's `NaifKeywords.BODY_RADII`.

Rationale: ISIS `cam2map` uses the DEM by default for its projection,
and we validated that using the same DEM brings csm2map's coverage
agreement with ISIS from 96.5% (ellipsoid) to 99.95% (DEM). The DEM
path matters because of how CSM's iterative `groundToImage()`
converges differently when the surface model is a sphere vs a
topography-aware elevation field — especially at the image edges
where a few meters of relief can shift a "marginal" pixel's
mapping by a pixel or two.

The DEM is sampled bilinearly at the coarse-grid lat/lon points
(a few hundred thousand per image), and the results are passed as
per-point surface radii to `groundToImage`. A rasterio handle is kept
open across the whole run so subsequent windowed reads hit the FS
cache cleanly.

Consequence: csm2map matches ISIS `cam2map`'s valid-pixel mask to
99.95%+ on CTX validation cubes with no flag tweaking.

### 3.4 Numerical precision: when float32, when float64

This is the decision most likely to come up in review. The short
answer: the **pixel-coordinate arithmetic runs in float32** (to
halve memory bandwidth), but the **SPICE/geometric math stays in
float64** (because that's what NAIF and CSM use internally). The
design is deliberate and we measured the impact.

**What is float32:**

- `coord_map.input_lines` and `coord_map.input_samples` — the
  per-output-pixel "where in the input cube to look" arrays. These
  are pure pixel coordinates.
- Bilinear interpolation weights and intermediate arrays in
  `_bilinear_upsample_pair`.
- Input cube DN values (as they arrive from rasterio / np.fromfile).
- Output cube DN values.
- DEM radii (stored float32 after the int16 → meters conversion).

**What stays float64:**

- pyproj coordinate transforms (map ↔ lat/lon) — pyproj's C core
  computes in double internally and always returns doubles.
- CSM `groundToImage` / `imageToGround` — C++ CSM plugin, always
  computes in double. We convert to float32 only when storing the
  result for downstream use.
- `scipy.ndimage.map_coordinates` — its C implementation always
  computes the bicubic interpolant in double regardless of input
  dtype. We pass float32 in and get float32 out, but the arithmetic
  inside scipy is full-precision.
- ALE ISD generation — ALE writes and `csmapi.Isd` parses JSON with
  double-precision numerics throughout.

**Why float32 is safe for pixel coordinates on line-scan and framing
planetary cameras:**

The IEEE 754 float32 mantissa is 24 bits, giving roughly 7 decimal
digits of precision. For a pixel coordinate of magnitude N the
precision is `N × 2⁻²³`. Working out the numbers:

| Camera / context | Max coord | float32 precision at that magnitude |
|------------------|----------:|------------------------------------:|
| JANUS (2000 × 2000) | 2000 | 0.00024 px (~120 nm on the detector) |
| CTX (2536 × 12288) | 12288 | 0.0015 px (~9 µm detector) |
| HiRISE NAC (20,000 × 100,000) | 100,000 | 0.012 px (12 µm detector equivalent) |
| HRSC line-scan (5184 × very-long) | 200,000 | 0.024 px |

Compare to:

- the intrinsic scatter of the CSM iterative solver (~10⁻⁶ px)
- the intrinsic scatter of the underlying SPICE CK interpolation
  (typically 10⁻³ to 10⁻² px for reconstructed kernels)
- the stated accuracy of published CSM validation papers (~0.1 px
  framing, ~0.5 px line-scan)

Even for HiRISE — the most demanding case in planetary imaging —
float32 pixel coordinates are **two orders of magnitude below**
the camera model's own error floor. There is no measurable loss.

**Measured impact on the MRO CTX J08 validation run**:

| | Original (float64 path) | Current (float32 path) |
|---|---:|---:|
| Both valid | 27,357,115 | 27,357,117 |
| ISIS-only | 14,804 | 14,802 |
| CSM-only | 3,401 | 3,399 |
| Mean (CSM − ISIS) | 0.000000 | 0.000000 |
| Median | −0.000008 | −0.000008 |
| Std | 0.001101 | 0.001101 |
| Max \|diff\| | 0.015722 | 0.015720 |
| \|diff\| < 0.01 | 100.00% | 100.00% |

Two pixels flipped between "valid" and "invalid" at the footprint
boundary — out of 46,404,923 — where a marginal bounds check
(`in_lines <= input_n_lines - 0.5`) falls on one side or the other
depending on whether the value is computed in float32 or float64.
That's a 4×10⁻⁸ effect and it sits below the CSM vs ISIS camera
model disagreement itself (which is the dominant boundary noise
floor at ~14 K pixels). All DN statistics are identical to six
decimal places.

**Where we would NOT recommend float32:**

- A downstream pipeline that consumes `coord_map` directly and needs
  to detect sub-0.001-pixel discrepancies (e.g. rigorous co-register
  analysis). In that case the correct answer is to expose a
  `precision="float64"` option that keeps coordinates in doubles
  through the whole pipeline. We don't currently expose this because
  no use case requires it, but adding it is a ~20-line change.
- Very-very-long-baseline line-scan cameras (Earth remote sensing
  has some with >10⁶ lines). None of the planetary cameras we target
  are in this regime.

### 3.5 Coarse grid + vectorized bilinear upsample, `step=32` default

Decision: we evaluate CSM `groundToImage` on a coarse grid of every
`step`-th output pixel (default `step=32`), then bilinearly upsample
the resulting `(line, sample)` coordinate map to the full output
resolution using a hand-vectorized numpy broadcast. This is
essentially what ISIS's `ProcessRubberSheet` does with its quad-tree
of bilinear patches, but as a uniform grid that is trivially
vectorizable.

Rationale:

- A per-pixel CSM evaluation of every output pixel (the `--dense`
  path) is O(10⁸) CSM calls for a typical CTX output. At ~240 K
  CSM calls per second that is ~7 minutes per image. Not acceptable.
- The ground-to-image mapping is smooth over the length scales of
  typical line-scan instruments on smooth bodies: over a 32-pixel
  output cell (~192 m on the ground at 6 m/px) the mapping varies
  by much less than a pixel. Bilinear approximation has negligible
  error.
- We measured the accuracy explicitly: at `step=16`, `step=32` and
  `step=64` on the J08 CTX reference cube, the DN comparison
  against ISIS `cam2map` is identical to 4 decimal places
  (std=0.001101 at all three steps, max \|diff\| within rounding
  noise, pixel counts within ±3 of 27.36 M total). `step=64` saves
  another ~200 ms over `step=32` but `step=32` is a safer default
  because it leaves headroom for higher-resolution cameras where
  the assumption of local smoothness is less safe.

Exposed as `--step`. Users with HiRISE-class resolution (25 cm/px)
or images of rugged terrain near the limb should drop to `--step 16`
or `--step 8`.

Consequence: ~182 K CSM calls become ~46 K CSM calls. The
`coord_transform` stage drops from ~3.4 s (scipy
`RegularGridInterpolator`) to ~0.63 s (vectorized bilinear with
threaded stripe execution at float32).

### 3.6 Per-pixel bounds check AFTER interpolation

Decision: we do not mask out-of-bounds coarse-grid points (points
where CSM would map to `sample < 0` or `sample > n_samples` in the
input cube) before the bilinear upsample. Instead, we let the coarse
grid interpolate naturally and apply a per-output-pixel bounds check
to the interpolated result.

Rationale: If we NaN the out-of-bounds coarse points, the bilinear
interpolator propagates NaN through any coarse cell with a NaN
corner, which produces a ragged edge that is ~1 pixel coarser than
the true camera footprint. Applying the bounds check AFTER
interpolation gives us a per-pixel-accurate footprint boundary even
though the coordinate map itself was computed on a coarse grid.

Consequence: csm2map's valid-pixel mask tracks the true camera
footprint at per-pixel resolution despite the coarse-grid
optimization. Small but visible at the footprint edges in a visual
comparison with a naive coarse-grid implementation.

### 3.7 Grid alignment and snap rule

Decision: when the input MAP PVL file contains explicit
`UpperLeftCornerX`, `UpperLeftCornerY` and `Samples`/`Lines`
keywords, we use those verbatim. When it doesn't, we build the grid
from the lat/lon ground range using a simple snap rule (floor the
minimum X and Y to a multiple of the pixel resolution, round the
extent to the nearest integer pixel count).

Rationale: ISIS `cam2map` uses a version-specific snap rule that we
could not precisely replicate from the dev-branch source (the actual
installed ISIS 9.0.0 behavior differs by ~1-2 pixels from what the
dev source would produce). Rather than chasing version-specific
idiosyncrasies, we accept whatever grid the user gives us. For
pixel-perfect comparison against an existing ISIS output, the
recommended workflow is to read the output cube's Mapping group and
pass it (with its explicit UpperLeftCornerX/Y) back in as the MAP
file for the csm2map run.

Consequence: csm2map can produce a pixel-identical grid to any ISIS
output when given that output's own Mapping group as input, without
reverse-engineering the ISIS snap algorithm.

### 3.8 Output container: GeoTIFF with ZSTD compression

Decision: the output is a tiled float32 GeoTIFF compressed with
ZSTD level 3, 256×256 tile size, `NUM_THREADS=ALL_CPUS`.

Rationale:

- **GeoTIFF over ISIS cube**: GeoTIFF is readable by every GIS tool
  (QGIS, ArcGIS, GDAL, rasterio, gdal2tiles, mapserver, leaflet
  tile servers...) without an ISIS installation. ISIS cube is a
  dead end for non-ISIS downstream pipelines.
- **Tiled**: gives cloud-friendly random access, efficient for
  future COG (Cloud-Optimized GeoTIFF) workflows.
- **ZSTD over LZW**: ZSTD at level 3 is ~8× faster than LZW at write
  time on our benchmark (0.13 s vs 1.01 s for a 100 MB float32
  image) AND produces a smaller file (81 MB vs 103 MB) because
  ZSTD's dictionary compression is a better fit for smooth image
  data than LZW's byte-level patterns.
- **Compatibility**: ZSTD was added to the TIFF spec via libtiff
  4.0.10 (November 2018) under compression tag 50000. GDAL 2.3+
  handles it transparently. Every GDAL-based tool (QGIS, rasterio,
  ArcGIS Pro 2.3+, ENVI recent versions, GRASS, R terra,
  scikit-image via GDAL) can read ZSTD GeoTIFFs. The compatibility
  risk is limited to pre-2018 toolchains, which is an acceptable
  trade.

If a user needs maximum compatibility (e.g. for distribution to an
unknown audience) they can fall back to LZW or DEFLATE via a
`--compression` flag. We have not yet exposed this flag; it's on
the roadmap.

## 4. Validation

Two CTX test cubes have been used as the acceptance benchmark:

- `J08_048038_1842_XN_04N287W.lev1.cub` — 2536 × 12288, Syrtis Major
- `F09_039335_1833_XI_03N284W.lev1.cub` — 5000 × 7168, different
  orbit and aspect ratio

For each, we run `csm2map` and ISIS 9.0.0 `cam2map` on the same cube
with the same MAP file and compare pixel-for-pixel in the overlap
region. Acceptance criteria:

1. **DN agreement**: 100% of overlapping pixels within 0.01 DN
2. **Coverage**: ≥99.9% of ISIS valid pixels also in csm2map
3. **No systematic bias**: mean \|DN diff\| < 10⁻⁵
4. **Std**: DN std < 0.002 (≈ half-pixel bicubic ringing noise)

Results on both cubes:

| | J08 | F09 |
|---|---:|---:|
| Both valid | 27,357,117 | 31,471,845 |
| Coverage overlap | 99.95% | 99.96% |
| Mean (CSM − ISIS) | 0.000000 | 0.000000 |
| Std | 0.001101 | 0.001192 |
| Max \|diff\| | 0.015720 | 0.011567 |
| \|diff\| < 0.001 | 68.20% | 66.09% |
| \|diff\| < 0.01 | 100.00% | 100.00% |

All four criteria pass on both cubes. The ~14-18 K pixels of
coverage mismatch are the structural CSM vs CSPICE camera noise
floor described in §3.1 — see `scripts/disagreement_analysis.md`
for the detailed investigation.

## 5. Performance

Single-image wall-clock on an Apple Silicon laptop (8 cores),
warm filesystem cache, MRO CTX J08 cube (125 MB input, 46 M output
pixels, DEM-aware projection):

| Stage | Time | % |
|---|---:|---:|
| load_camera (ALE → ISD → RasterGM) | 0.62 s | 22% |
| target_radii | 0.10 s | 4% |
| dem_open (windowed lazy) | 0.12 s | 4% |
| build_grid | 0.01 s | <1% |
| coord_transform (coarse CSM + bilinear upsample) | 0.63 s | 23% |
| read_input (cube → numpy) | 0.15 s | 5% |
| resample (threaded scipy map_coordinates) | 1.00 s | 36% |
| write_output (ZSTD GeoTIFF) | 0.13 s | 5% |
| **Total** | **~2.8 s** | **100%** |

Comparison to ISIS `cam2map` (sandbox, `ulimit -n 4096` for the
overflow fix documented in `scripts/isis_sandbox_fix.md`): ~65 s.
Speedup: **~23×**.

Performance headroom not yet used:

- ALE ISD caching with mtime-based invalidation (−0.6 s per run on
  warm cache, ~22% of wall time)
- Overlapping `read_input` with `coord_transform` via a background
  thread (−0.15 s, ~6%)
- GPU port of the resample and/or the bilinear upsample (−0.5 to
  −1.0 s, 20-35%) — deferred until the CPU path hits a real wall

## 6. Known limitations

1. **CSM vs CSPICE camera boundary noise**: ~0.04% of pixels at the
   footprint boundary disagree between csm2map and ISIS `cam2map`.
   Irreducible without using ISIS's exact camera class. Documented.
2. **No bundle adjustment**: csm2map projects with whatever pointing
   is in the cube's embedded blobs. If those are pre-jigsaw, the
   projection uses pre-jigsaw geometry. The `spice_source=isis`
   default makes this correct after jigsaw, but csm2map itself does
   not do BA.
3. **Rolling-shutter-aware CSM models are supported, jitter
   coefficients are not yet exposed**: `UsgsAstroLineScanSensorModel`
   can consume jitter coefficients in its ISD, but we currently
   don't expose any knobs to override them. For framing cameras with
   rolling shutter (Europa Clipper EIS NAC, JANUS) the default
   behavior is correct as long as the ALE driver produces the right
   ISD.
4. **Multi-band handling is shallow**: the current code treats
   multi-band cubes as independent single-band projections. No
   band-to-band geometric registration tricks are applied. Fine for
   CTX (single band), probably fine for JANUS filter-wheel framing
   mode.
5. **No stereo / DEM generation**: `csm2map` is a projection tool
   only. Stereo DEMs still require ASP.
6. **No Planetocentric/Planetographic conversion when the MAP file
   disagrees with the CSM model**: ISIS will warn you and convert;
   we currently trust the MAP file. This hasn't bitten us on CTX
   because the MAP files are always planetocentric but it's a
   latent bug waiting for a planetographic MAP input.
7. **`--clip-to-footprint` does not match ISIS `cam2map`** (historical
   note — the flag was added under a hypothesis that turned out to
   be wrong). Early in development we saw a ~3% coverage gap vs
   `cam2map` and observed that ~99.7% of the excess pixels sat
   outside the `footprintinit` polygon stored in the cube. We added
   `--clip-to-footprint` on the working assumption that `cam2map`
   was internally clipping to that polygon. A later direct source
   read of ISIS 9.0.0 `cam2map.cpp` showed no reference to any
   polygon at all, and an empirical test (stripping the polygon
   from the cube and re-running `cam2map`) produced bit-identical
   output to the unmodified run. The original coverage gap was
   actually an ellipsoid-vs-DEM shape model difference; once we
   integrated MOLA as the default shape model, the gap closed to
   99.95% **without any polygon clipping**. The flag is retained in
   0.7.x only as an escape hatch for downstream tooling that wants
   a polygon-shaped output mask; it is explicitly *not* a
   compatibility mode for `cam2map`. See `docs/csm2map.qmd §
   "The footprintinit polygon precision story"` for the full
   narrative; it's worth including in the paper as a concrete
   example of a hypothesis-test-reject debugging loop during
   validation.

## 7. Non-goals

- Replacing ISIS. csm2map replaces `cam2map`, not the whole ISIS
  ecosystem. `jigsaw`, `qmos`, `findimageoverlaps`, `autoseed`, etc.
  continue to live in ISIS and that's fine. If someone wants an
  ISIS-free pipeline end-to-end they need the hybrid-pattern
  prototype (see `Plans/hybrid-pattern-prototype.md`), not just
  csm2map.
- Being a drop-in binary replacement for `cam2map`. We match the
  behavior that matters (DN values, coverage, camera geometry) but
  we produce GeoTIFF, not ISIS .cub, and we don't implement every
  `cam2map` parameter.
- Generalizing to arbitrary CSM models we haven't tested. We only
  claim validity for the CSM model classes usgscsm exposes
  (frame, line scan, push frame, SAR), and we only claim sub-pixel
  correctness on the instruments we've explicitly validated (MRO
  CTX). Extension to other instruments requires the
  instrument-specific validation dance documented in the
  companion `csm_research.md` note.

## 8. Future work

Tracked as separate plan documents in `Plans/`:

- `Plans/hybrid-pattern-prototype.md` — controlled multi-image
  mosaicking via arosics + a CSM-native pose refiner, ISIS-free
- `Plans/gpu-acceleration.md` — torch MPS / cupy backends for
  resample and bilinear upsample, if the CPU path ever becomes the
  bottleneck

Smaller items worth doing in the near term:

- `--compression {zstd,deflate,lzw,none}` CLI flag for users who
  need legacy-tool compatibility
- `--precision {float32,float64}` CLI flag for sub-pixel-accuracy
  downstream use cases that we don't have today but might later
- ALE ISD caching with mtime-based invalidation (saves ~0.6 s per
  repeated run on the same cube)
- JANUS / Europa Clipper EIS validation once real (or simulated)
  data becomes available

## 9. Citation

When citing `csm2map` in publications, please cite:

- The `isistools` package DOI (TBD — will be registered at the
  time of paper release)
- This design document (reference the commit SHA it was read at, or
  the Zenodo archive if released separately)
- The CSM and ALE papers that we build on:
  - Laura, Mapel & Hare 2020, *Earth and Space Science* 7(6),
    doi:10.1029/2019EA000713
  - Laura et al. 2024, *Remote Sensing* 16(4):648,
    doi:10.3390/rs16040648

## 10. References

- NGA CSM TRD v3.0: <https://gwg.nga.mil/documents/csmwg/documents/CSM_TRD_Version_3.0__15_November_2010.pdf>
- usgscsm: <https://github.com/DOI-USGS/usgscsm>
- ALE: <https://github.com/DOI-USGS/ale>
- swigcsm (Python bindings): <https://github.com/DOI-USGS/swigcsm>
- ISIS3: <https://github.com/DOI-USGS/ISIS3>
- ESA SPICE for JUICE (JANUS): <https://www.cosmos.esa.int/web/spice/spice-for-juice>
