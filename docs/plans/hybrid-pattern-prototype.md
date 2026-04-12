# Plan: Hybrid CSM + arosics controlled-mosaicking prototype

## Context

We want a controlled multi-image mosaicking pipeline that is fully
ISIS-free and CSM-native, as the reference implementation for JANUS
(and any future mission where ISIS doesn't exist yet). The pipeline
combines:

- `csm2map` (ours) for physically-correct initial map projection using
  CSM camera models
- `arosics` (pip package, GFZ Potsdam) as a robust sub-pixel tie-point
  detector against a reference basemap
- A lightweight CSM-native pose refiner (scipy.optimize, not ASP
  `bundle_adjust`, so no ISIS dependency at all)
- `csm2map` again with the refined CSM states for the final projection

The architectural bet: arosics is excellent at sub-pixel image matching
but does 2D rubber-sheet warping that isn't physically correct. A
bundle adjuster is physically correct but heavy. By using arosics as a
*tie-point detector* feeding a *small bundle adjuster*, we get the
sub-pixel matching accuracy of arosics with the physical-model rigor
of bundle adjustment, and no ISIS.

This is explicitly a prototype. The goal is to prove the pattern works
end-to-end on a real CTX test case before investing in a production
implementation for JANUS (where no data exists yet anyway).

## Scope (explicit non-goals)

- **Not building a general-purpose bundle adjuster.** The pose refiner
  in this prototype handles framing cameras and simple line-scan cases
  only, just enough to demonstrate the pattern on CTX. A full BA with
  jitter solving, rolling-shutter handling, and sparse constraint
  systems is a production project, not a prototype.
- **Not replacing `csm2map`.** This is a pipeline that composes
  `csm2map` with other tools; `csm2map` itself stays as it is.
- **Not adding ISIS anywhere.** The whole point is to prove the
  ISIS-free path. If we hit something that needs ISIS, we either work
  around it or mark it as a blocker and scope it out.
- **Not a user-facing CLI on day one.** Prototype is a Python script
  or Jupyter notebook. CLI comes later if the pattern works.

## Pipeline

```
  ┌─────────────┐
  │ N CTX cubes │
  └──────┬──────┘
         │
         │ csm2map (current code, unchanged)
         │  reads cube blobs → CSM ISD → usgscsm model → GeoTIFF
         ▼
  ┌─────────────┐
  │ N rough     │      Projected in map coords using initial pointing
  │ projected   │      from the cube's embedded SPICE (or jigsaw-refined
  │ GeoTIFFs    │      if you ran jigsaw earlier, thanks to spice_source=isis).
  └──────┬──────┘
         │
         │ arosics COREG_LOCAL vs reference basemap
         │  - basemap: Caltech CTX Global Mosaic for CTX validation
         │             or any controlled basemap at ≥ target resolution
         │  - output: grid of tie points per image in (input_px, lat, lon)
         ▼
  ┌─────────────┐
  │ tie points  │      Per-image list of (image_sample, image_line,
  │ per image   │      ref_lat, ref_lon, residual)
  └──────┬──────┘
         │
         │ scipy-based CSM pose refiner
         │  - cost: image-space residual at each tie point
         │  - variable: 6-DOF pose correction per image
         │    (or just attitude, depending on CSM state schema)
         │  - solver: least_squares with trust region
         │  - per image or joint (all images, shared tie points) depending on scale
         ▼
  ┌─────────────┐
  │ refined     │      Updated CSM states, saved as JSON next to cubes
  │ CSM states  │
  └──────┬──────┘
         │
         │ csm2map with refined states
         │  - needs a new flag: --csm-state to bypass ALE and use a
         │    pre-computed CSM state JSON directly
         ▼
  ┌─────────────┐
  │ N controlled│      Now registered to the basemap frame,
  │ GeoTIFFs    │      suitable for mosaicking with gdalbuildvrt +
  └─────────────┘      gdalwarp or rasterio.merge
```

## Components to build

### 1. CSM state extraction and injection in csm2map

Needed because the pose refiner has to be able to (a) read the current
CSM state of a model so it has a starting point to optimize, and (b)
write a refined state back that `csm2map` will use on the next run.

Current `csm2map/processing/camera.py::load_camera()` builds the CSM
model via ALE → ISD JSON → `plugin.constructModelFromISD()`. We need a
second path:

```python
def load_camera_from_state(state_json_path: Path) -> csmapi.RasterGM:
    """Load a CSM model from a pre-computed state JSON, bypassing ALE."""
    ...
```

using `plugin.constructModelFromState()` instead of `constructModelFromISD`.
usgscsm supports this natively — it's how refined models get persisted.

Also add a state save helper:

```python
def save_camera_state(model: csmapi.RasterGM, out_path: Path) -> None:
    """Dump the CSM model state JSON so a refined model can be reloaded."""
    state = model.getModelState()  # returns a JSON string
    out_path.write_text(state)
```

And plumb a `--csm-state PATH` option through the CLI so
`isistools csm2map --csm-state refined.json input.cub out.tif` uses
the refined model instead of regenerating one from the cube.

### 2. arosics tie-point extractor

A small Python function that takes a projected image (GeoTIFF) and a
reference basemap (GeoTIFF), runs `arosics.COREG_LOCAL`, and returns
the tie-point grid as (image_sample, image_line, ref_lon, ref_lat)
tuples.

```python
def extract_tie_points(
    image_geotiff: Path,
    basemap_geotiff: Path,
    grid_step: int = 50,
    max_shift_px: float = 5.0,
    min_reliability: float = 70.0,
) -> np.ndarray:  # shape (N, 5): [samp, line, lon, lat, residual]
    ...
```

Implementation notes:

- arosics `COREG_LOCAL` naturally produces this grid internally; we
  just extract it rather than applying the warp
- use `max_shift_px` to reject obvious mismatches
- use `min_reliability` (arosics correlation coefficient threshold)
  to filter low-confidence points
- convert tie points from *output projected-image coordinates* back
  to *input camera (sample, line) coordinates* using the same CSM
  model that produced the projection. This is the inverse of what
  csm2map does and needs a helper function

### 3. scipy pose refiner

For each image, optimize the CSM state to minimize image-space
residuals at the tie points:

```python
def refine_pose(
    model: csmapi.RasterGM,
    tie_points: np.ndarray,  # (N, 5): image_samp, image_line, lon, lat, _
    mode: str = "attitude_only",  # or "6dof"
) -> csmapi.RasterGM:
    """Return a refined model whose groundToImage at the reference lat/lon
    best matches the observed image_sample/image_line at each tie point."""
    ...
```

Cost function: for each tie point, convert (lon, lat) → ECEF → call
`model.groundToImage()` → compare to observed (image_samp, image_line)
→ accumulate residual. Use `scipy.optimize.least_squares` with the
trust region reflective method.

Variable dimensions:

- `attitude_only` mode: 3 params per image (small Euler perturbation
  of the pointing quaternion). This is what ISIS jigsaw's default is,
  and what most planetary missions actually need.
- `6dof` mode: 6 params per image (3 attitude + 3 position). More
  flexible, more underdetermined, needs more tie points.
- For the prototype: start with `attitude_only`, add `6dof` later if
  needed.

Per-image refinement first. Joint refinement across N images with
shared ground points is a second-phase improvement and much harder.
Explicit non-goal for v1.

### 4. Glue script

A top-level Python script (not yet a CLI) that drives the pipeline:

```python
# hybrid_controlled_mosaic.py
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("cubes", nargs="+", type=Path)
parser.add_argument("--basemap", type=Path, required=True)
parser.add_argument("--map", type=Path, required=True, help="MAP PVL file")
parser.add_argument("--outdir", type=Path, required=True)
args = parser.parse_args()

for cube in args.cubes:
    # 1. Initial projection
    rough_tif = args.outdir / f"{cube.stem}_rough.tif"
    run_csm2map(cube, rough_tif, map_file=args.map)

    # 2. Tie points vs basemap
    tie_points = extract_tie_points(rough_tif, args.basemap)

    # 3. Pose refinement
    model = load_camera(cube, spice_source="isis")
    refined = refine_pose(model, tie_points)

    # 4. Save refined state
    state_path = args.outdir / f"{cube.stem}.csm_state.json"
    save_camera_state(refined, state_path)

    # 5. Re-project with refined pose
    final_tif = args.outdir / f"{cube.stem}_controlled.tif"
    run_csm2map(cube, final_tif, map_file=args.map, csm_state=state_path)

# 6. Mosaic
run_gdalbuildvrt([...], args.outdir / "mosaic.vrt")
run_gdalwarp(args.outdir / "mosaic.vrt", args.outdir / "mosaic.tif")
```

## Validation test

Use the existing J08 and F09 CTX cubes plus a couple more adjacent CTX
frames to build a small 3–4 image Syrtis Major mosaic. Compare:

- Baseline: direct `csm2map` projections stitched without refinement
- Arosics-direct: arosics COREG_LOCAL applied after csm2map,
  warped outputs stitched (no bundle adjustment — pure 2D warp)
- Hybrid: the pipeline above, with refined CSM states feeding csm2map

Metric: overlap residuals at seam lines. For adjacent frames, compute
the DN difference in the overlap zone as a proxy for geometric
accuracy. A well-controlled mosaic should have near-zero seam
residuals in low-relief areas.

Reference basemap: Caltech/JPL CTX Global Mosaic (Dickson et al.) at
~5 m/pix. Download the relevant tiles for the Syrtis Major region.

## Dependencies

- `arosics` (new, pip install)
- `scipy.optimize.least_squares` (already have)
- `rasterio` (already have)
- `numpy` (already have)
- `usgscsm.getModelState()` / `constructModelFromState()` — need to
  verify these are exposed in the Python bindings we built

No new conda packages, no new C++ builds, no ISIS dependency.

## Open questions to resolve during implementation

1. Does the Caltech CTX Global Mosaic come tiled at a usable size for
   arosics's windowed correlation, or do we need to crop to the
   relevant region first?
2. What's the exact schema of `UsgsAstroLineScanSensorModel`'s model
   state JSON? We need to identify which fields correspond to
   attitude/position so the refiner knows what to optimize.
3. Can `constructModelFromState` round-trip a state that was produced
   by `getModelState`, or do we need to massage the JSON? (ISIS Jigsaw
   does some massaging — might need the same.)
4. For attitude-only refinement, do we modify the first quaternion
   (image start) or all per-line quaternions? For framing cameras it's
   one quaternion; for line-scan it's a function of time and the
   refinement is more subtle.
5. How many tie points does arosics typically produce on flat Mars
   terrain at 6 m/pix? Need enough to over-constrain a 3-DOF
   attitude fit (so ≥6, practically ≥30 for a decent least-squares).

## Verification steps

- [ ] `load_camera_from_state` / `save_camera_state` round-trip
  produces a bit-identical model (run `csm2map` with both and compare
  output)
- [ ] arosics tie points on two adjacent CTX frames, visualized,
  look reasonable (not scattered random, clustered in high-texture
  regions)
- [ ] Pose refinement on a single cube against itself (zero-residual
  case) returns unchanged pose
- [ ] Pose refinement on a cube with a known injected attitude offset
  recovers the offset to within 0.1 px
- [ ] Mosaic seam residuals drop measurably between baseline and
  hybrid modes

## Estimated effort

~2–3 days of focused work for the prototype, assuming no surprises
with usgscsm state serialization. Most of the code is straightforward
glue; the pose refiner is the only piece with real numerical content.
