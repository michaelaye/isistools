# Plan: GPU acceleration for csm2map

## Revision history

- **2026-04-11** — first draft, when csm2map was a ~11 s/image prototype
  on J08, proposing "thread the CSM loop" as the biggest quick win.
- **2026-04-12** — revised after the 0.8.0 body-agnostic refactor and a
  fresh F05 re-profile. The "thread the CSM loop" quick win is already
  shipped (see §4 below) and the baseline numbers from the first draft
  are stale. GPU resample is still the single biggest improvement
  available.

## Context

By v0.8.0 csm2map runs a 1 GB full-length F05 CTX cube through the full
pipeline (camera → coord transform → resample → ZSTD GeoTIFF) in
**~44.6 s** warm cache, producing a 557M-pixel output. That's a ~13×
speedup over ISIS `cam2map` on the same inputs (556 s on the same
machine, measured in the same session).

For one-off work 44 s is fine. For bigger campaigns or for JANUS /
HiRISE scale it's worth a further 2–3× speedup. This plan identifies
where the time actually goes **today** (not in the stale 0.7.0
numbers the first draft was written from), what the realistic ceiling
is, and how much engineering it costs to get there.

## Current profile (0.8.0, F05 warm cache, `--profile`)

Averaged over three runs, warm cache (first cold run is ~8 s slower in
`coord_transform` due to DEM windowed-read paging — discarded):

| Stage | Wall time (s) | % of total | GPU-addressable? |
|-------|---------------|------------|------------------|
| ALE ISD + camera load | 1.28 | 2.9% | No — Python + SPICE |
| DEM open (window cache warm) | 0.13 | 0.3% | No |
| Grid build | 0.01 | 0.0% | No |
| **`coord_transform`** (coarse CSM + bilinear upsample) | **12.2** | **27.3%** | **Partial** — only the upsample, not the CSM calls |
| Read input cube | 0.71 | 1.6% | No — IO bound |
| **`resample`** (scipy `map_coordinates`) | **28.7** | **64.3%** | **Yes — biggest win** |
| GeoTIFF write (ZSTD, threaded) | 1.50 | 3.4% | No — IO + GDAL |
| (other) | 0.11 | 0.3% | — |
| **total** | **44.64** | **100%** | |

**Two things stand out vs the 0.7.0 draft:**

1. **Resample is now 64% of wall time**, not 27%. The float32 + threaded
   resample path in 0.7.0 made the stage faster in absolute terms but
   everything else got faster with it, so resample's *share* grew.
   GPU-accelerating resample has a bigger absolute payoff now than it
   did in the 0.7.0 draft.

2. **`coord_transform` is a single stage of 12.2 s** that combines the
   coarse CSM calls (~10 s), the DEM sampling (~0.1 s), and the
   vectorized float32 bilinear upsample (~2 s). The upsample is the
   only GPU-addressable part; the CSM calls themselves are a C++
   sensor-model root-finder that can't run on GPU without an upstream
   port of `usgscsm`.

## What moved since the 0.7.0 draft (and what was wrong)

1. **"Quick win: thread the CSM ground_to_image loop"** was proposed in
   the 0.7.0 draft as a trivial `ThreadPoolExecutor` change. That code
   is now in `processing/camera.py::ground_to_image_batch` with
   `workers = os.cpu_count()` as the default. **In practice it gives
   roughly zero speedup** because csmapi's SWIG bindings hold the GIL
   inside `groundToImage`. The scaffolding is retained so that if a
   future csmapi release marks `groundToImage` as nogil (the standard
   way C++ extensions release the GIL), the speedup materializes
   automatically. But **do not expect threading alone to move the
   needle on the CSM call stage**. The 0.7.0 draft's claim of 4–8×
   from this quick win was wrong.

2. **Float32 throughout the pipeline** was listed as a future item in
   the draft. It's now done and already contributes to the baseline
   above. No more free speedup to capture there.

3. **`_bilinear_upsample_pair` vectorized + stripe-threaded** in
   `processing/transform.py` — also done, also in the baseline. The
   previous `scipy.ndimage.zoom` call was the bottleneck of the upsample
   sub-stage and is now gone.

4. **ZSTD with `num_threads=ALL_CPUS`** for the writer. Already the
   fastest portable option for big tiled GeoTIFFs; not a GPU target.

5. **ISD JSON caching (`--cache-isd` flag)** is still not implemented.
   It would save the ~1.28 s `load_camera` stage on repeated runs over
   the same cube (handy for benchmarks and dev loops) but it's only
   ~2.9% of current wall time — marginal, noted, low priority.

## The one GPU win that's actually worth pursuing: resample

scipy `map_coordinates` at order=3 (bicubic) on a 557M float32 output
is memory-bandwidth-limited and single-GPU-kernel territory.

### Option A — `torch.nn.functional.grid_sample` on MPS (macOS) / CUDA (Linux)

Single backend, runs everywhere the maintainer does. The trade-off:
torch's `grid_sample` is NOT the same API as scipy's `map_coordinates`.

Key differences to reconcile:

- **Coordinate convention**: `grid_sample` uses normalized coordinates
  in `[-1, 1]` referring to the input image corners, not raw pixel
  indices. Our `coord_map` currently stores raw (line, sample) pixel
  indices. Conversion: `norm_x = 2 * sample / (n_samples - 1) - 1`,
  same for `norm_y`. One-time helper.
- **Out-of-bounds behavior**: `grid_sample` defaults to `padding_mode=
  "zeros"`. Our code uses scipy's `mode="constant", cval=fill_value` —
  which is effectively the same when `fill_value=0`, different when
  we later set `fill_value=np.nan`. We apply the mask *after* resample
  anyway (`result[~coord_map.valid] = fill_value`), so the torch
  default is fine.
- **Interpolation modes**: torch supports `bilinear` and `bicubic` —
  matches our `Interpolation.BILINEAR` and `Interpolation.BICUBIC`.
  For `NEAREST`, torch has `mode="nearest"`. Full coverage.
- **Precision**: torch defaults to float32, which matches our existing
  float32 intermediate. MPS on Apple Silicon has solid float32
  support; float64 is slow/unsupported on MPS. We already float32.

**Expected performance**: roughly 10× over scipy's threaded CPU path on
a memory-bandwidth-bound kernel. On F05's 557M pixels that would take
the stage from 28.7 s to ~3 s, shaving ~26 s off total pipeline time.

**End-to-end speedup estimate**:
- Current: 44.6 s
- With GPU resample only: ~18.9 s → **~2.4× overall**
- Ceiling: limited by `coord_transform` (12.2 s, mostly CSM calls that
  can't move to GPU without upstream work).

### Option B — `cupy.ndimage.map_coordinates`

1:1 API match to scipy (drop-in replacement), but **cupy has no
Apple Silicon support** and never will unless AMD/NVIDIA ship on ARM
Mac hardware. This forces a second backend on Linux only, duplicating
the resample module for a dev environment that's specifically not
the maintainer's primary platform.

**Decision**: skip cupy. Go straight to torch MPS as the single GPU
backend — it runs everywhere, the coordinate-convention reconcilation
is a one-shot cost, and the rest of the code (coord map generation,
mask application) stays numpy.

### Option C — `jax.scipy.ndimage.map_coordinates`

Similar API to scipy, runs on CUDA + macOS Metal. Metal backend is
experimental. We'd inherit the JAX tracing/JIT machinery for one
function, which is overkill. **Skip** — torch is mature, JAX Metal is
not, and we don't need autodiff.

## What's NOT worth GPU-porting

1. **ALE ISD generation** — Python + SPICE kernel reads. Not a GPU
   workload. 2.9% of wall time, too small to matter.
2. **CSM `groundToImage` root-finder** — C++ in usgscsm, single-threaded
   per call, iterative Newton solver over a spacecraft trajectory.
   Needs a full CUDA port of `UsgsAstroLineScanSensorModel` to move,
   which is a person-month project that belongs in upstream usgscsm,
   not here. The upstream project is CPU-only and shows no sign of a
   GPU port. Note this as a longer-term ceiling; don't attempt it in
   this plan.
3. **DEM windowed reads** — rasterio/GDAL, dominated by decompression
   and memory allocation. Not GPU-addressable without refactoring the
   MOLA DEM itself (2 GB cube stored as int16) to a GPU-friendly
   tiled texture format. Out of scope.
4. **GeoTIFF write** — GDAL ZSTD, already multi-threaded. Compression
   itself is CPU work; GPU ZSTD exists but is niche and not plumbed
   into rasterio.
5. **`coord_transform`'s bilinear upsample** — already vectorized in
   numpy at float32 with thread striping. Moving it to GPU saves at
   most ~2 s and entangles the rest of the stage (which is the
   un-moveable CSM call) with GPU memory management. Don't bother.

## Implementation strategy (revised 2026-04-12)

Do these in order. Stop when the speedup is good enough for your use
case — you may not need the later steps.

1. **Prove-of-concept torch MPS resample** on a synthetic 10000×10000
   pixel test. Verify:
   - Sub-pixel accuracy matches scipy `map_coordinates` at order=3
     to better than 1e-4 in DN.
   - Coordinate-convention conversion (raw pixel → normalized) is
     reversible and correct.
   - float32 round-trip (numpy → torch → numpy) costs no more than
     ~0.3 s per stage transition on F05-scale arrays (~2.1 GB).
   - Effort: half a day.

2. **Backend abstraction in `processing/resample.py`**. Introduce a
   `_resample_band_gpu` function behind a `try/except ImportError`
   for `torch`. Default dispatch: GPU if available, CPU fallback
   otherwise. The `resample()` public API stays unchanged, and the
   package remains fully functional on machines without torch.
   Effort: half a day.

3. **End-to-end F05 benchmark**. Run csm2map --profile on F05 with
   the torch MPS path and compare against the 44.6 s baseline above.
   Target: **≤ 20 s total** (2.2×+ overall speedup). If the GPU
   stage comes in at ~3 s but total wall time sits above 25 s,
   there's a host↔device copy cost we didn't anticipate — investigate
   and tune before declaring success.
   Effort: half a day.

4. **Cold-cache instrumentation**. Measure the cold-vs-warm gap on
   `coord_transform` in the new pipeline. The current gap (8 s on
   cold, 0 on warm) is mostly DEM page cache. See if preloading the
   DEM bounding box on startup closes the gap — it's a free win
   unrelated to GPU work.
   Effort: half a day.

5. **ISD caching** (`--cache-isd` flag). Only a 1.3 s win but it's a
   1-function change in `camera.py::load_camera` and it lets
   benchmark scripts run the projection path in isolation. Effort:
   half a day.

**Total optimistic effort**: 2–3 days to ship a torch MPS backend
that halves end-to-end wall time on F05-class inputs, plus minor
cleanups. **Total pessimistic effort**: 4–5 days if the torch
coordinate conversion has subtleties we don't foresee.

## Expected outcome (0.8.0 baseline)

Best case (torch MPS resample + cold-cache fix + ISD caching, second
run onwards):

| Stage | Before | After | Notes |
|-------|-------:|------:|-------|
| load_camera | 1.28 | 0.15 | ISD cache hit |
| dem_open | 0.13 | 0.05 | preloaded |
| build_grid | 0.01 | 0.01 | unchanged |
| coord_transform | 12.2 | 10.0 | upsample on GPU saves ~2 s |
| read_input | 0.71 | 0.71 | unchanged |
| **resample** | **28.7** | **~3.0** | torch MPS `grid_sample` |
| write_output | 1.50 | 1.50 | unchanged |
| **total** | **44.6** | **~15.4** | **~2.9× overall** |

That's the realistic ceiling without touching the CSM call path.
Anything beyond that requires upstream GPU work in usgscsm itself.

## Non-goals

- **Not porting CSM / usgscsm to GPU.** Upstream project, wrong scope.
- **Not maintaining a cupy backend** alongside torch. Apple Silicon is
  the maintainer's primary platform; cupy doesn't run there.
- **Not training neural-net surrogates** for the sensor model. Cute,
  overengineered, accuracy unknown.
- **Not adding autodiff.** Nothing in csm2map needs gradients today.
- **Not requiring a GPU.** All GPU code paths must fall back to the
  existing CPU path, and the package must remain fully functional on
  CPU-only systems. This is a hard constraint — enforced by the
  `try/except ImportError` guard in step 2 above.

## Open questions

1. **How often does csm2map actually run per day?** If it's 5–10 times,
   a 44 s → 15 s speedup saves ~5 minutes per day — nice but not
   transformative. If csm2map is the inner loop of a nightly
   mosaicking job across 200 CTX cubes (~2.5 hours → ~50 minutes),
   this changes the feel of the tool. The answer decides whether the
   3-day torch port is worth it or if the current 0.8.0 CPU baseline
   is already good enough.

2. **Does the JANUS use case actually need GPU resample?** JANUS is a
   framing camera (2000×2000 × 13 filters). Output pixel count per
   image is ~50× smaller than F05. The CPU baseline may already be
   sub-second for JANUS. **Profile JANUS first before starting GPU
   work** — it might be a solved problem on CPU.

3. **MPS float32 stability on large tensors**. We'd be pushing ~2.1 GB
   of float32 through `grid_sample` in one shot. Apple Silicon unified
   memory handles this, but MPS has historically had quirks around
   large-tensor precision. The proof-of-concept in step 1 should test
   this directly.

4. **Is the 2–3× headline speedup actually worth 3 days of engineering
   time vs other things on the roadmap** (hybrid arosics+CSM
   prototype, JANUS enablement, pose-refinement hooks)? Compare
   against the `Plans/hybrid-pattern-prototype.md` plan — that one is
   blocked on CSM state I/O (1 day's work) and unlocks a new
   capability, while this plan unlocks ~2× on an existing one.

## Prioritization recommendation (for the roadmap, as of 2026-04-12)

If the question is "what's the single biggest ROI thing to do next":

1. **Hybrid prototype prerequisites** (CSM state load/save in
   `camera.py`) — 1 day, **unlocks a new capability** (controlled
   mosaicking). Higher ROI per engineering hour than GPU resample
   because capability expansion > 2× speedup on existing code.
2. **JANUS-first profile** — half a day, decides whether to do GPU
   work at all. If JANUS is already sub-second on CPU, GPU resample
   is a CTX-only optimization and its ROI drops.
3. **Then GPU resample** — 2–3 days, 2–3× overall speedup on CTX F05.

This plan is kept up to date as step 3 in the sequence. Revisit
this section after steps 1 and 2 actually happen.
