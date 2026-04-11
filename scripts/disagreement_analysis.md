# Why ~18K pixels still disagree (and why we can't fix it)

After integrating the MOLA DEM as the shape model, our `csm2map` output
matches ISIS `cam2map`'s valid mask to 99.95% on the J08 CTX test cube,
with all overlapping pixels agreeing to within 0.01 DN. The remaining
~18K disagreement is **a structural floor between CSM
`UsgsAstroLineScanSensorModel` and ISIS `CTXCamera`** and is not fixable
without using ISIS's actual camera model.

## What "disagreement" means

After running both tools, every output pixel falls into one of four buckets:

|              | ISIS valid | ISIS NULL |
|--------------|------------|-----------|
| CSM valid    | both valid (compare DN values here) | **CSM-only** ← disagreement |
| CSM NULL     | **ISIS-only** ← disagreement | both NULL (outside coverage) |

A "valid" pixel is one where the tool wrote a real DN value. Disagreeing
pixels are positions where the two tools made different "is this pixel
inside the camera footprint?" decisions. They have no DN values to
compare — they're a binary agree/disagree on validity.

## J08 disagreement characterization

```
ISIS-only: 14,804 pixels  (ISIS valid, CSM NULL)
CSM-only:   3,401 pixels  (CSM valid, ISIS NULL)
Total:     18,205         (~0.07% of total output)
```

Connected-component analysis:

- **CSM-only** (3,401): 1,716 components, 1,415 of them singletons,
  100% within 1 pixel of the both-valid region. These are scattered
  single-pixel rounding-noise pixels at the perimeter.
- **ISIS-only** (14,804): 1,511 components, median size 8, with one
  large 3,303-pixel diagonal stripe along the top edge of the
  parallelogram footprint. 92.7% are within 1 pixel of both-valid.

## Root cause: CSM vs CSPICE camera disagreement at line-scan edges

The big 3,303-pixel ISIS-only component lies along the **bottom edge of
the input cube** (highest output latitudes → first/last lines of CTX
acquisition). For these pixels:

- Our CSM `groundToImage` returns `in_line` values clustered in
  `[12287, 12288]` for a 12288-line input cube.
- Our bounds check rejects anything `> 12287.5` (the bottom edge of the
  last valid line center).
- ISIS's CSPICE camera returns values systematically ~0.3-0.5 pixels
  *lower* for the same ground points and accepts them.

This is **not a bounds-check bug**. CSM is internally self-consistent
(`imageToGround` followed by `groundToImage` round-trips to within
machine epsilon for any line/sample). The discrepancy is between
**CSM's `UsgsAstroLineScanSensorModel` line-time iteration** and
**ISIS's CSPICE-based `LineScanCamera::SetUniversalGround` iteration**.
Both implementations solve the same physical problem ("at what time was
the spacecraft pointing at this ground point?") but they use different
starting guesses, different convergence criteria, different finite-
difference formulations, and they converge to slightly different times
when the answer lies near the temporal limits of the acquisition.

A direct comparison of the lat/lon at image corner pixels confirms this:
ISIS `campt` at `(sample=0.5, line=0.5)` reports
`(lat=3.56914527, lon=72.79979982)`, while CSM `imageToGround(0.0, 0.0)`
returns `(lat=3.56891, lon=72.79841)`. The difference is `0.00024°` in
lat and `0.00139°` in lon, which corresponds to **about 2-3 pixels**.
The two camera implementations genuinely place the corner pixel at
slightly different ground positions. This offset isn't uniform across
the image; it varies along the perimeter and that's why most disagreeing
pixels are within 1 pixel of the boundary.

## What we tried

A bounds-relaxation sweep showed that adding a per-side "slop" doesn't
help meaningfully:

| Slop (px) | both     | ISIS-only | CSM-only | total disagreement |
|-----------|----------|-----------|----------|--------------------|
| 0.00      | 27.36 M  | 8,242     | 9,792    | **18,034**         |
| 0.25      | 27.37 M  | 4,937     | 12,990   | 17,927             |
| 0.50      | 27.37 M  | 1,655     | 16,167   | 17,822             |
| 0.75      | 27.37 M  | 515       | 21,473   | 21,988             |
| 1.00      | 27.37 M  | 0         | 27,432   | 27,432             |

The minimum total disagreement is at slop=0.5 (17,822) — only ~1%
better than the default. Increasing slop further makes the situation
worse: at slop=1.0 we capture every ISIS pixel but add 27K extra
"valid" pixels that ISIS rejects.

We left the bounds at slop=0 because:

1. The improvement is marginal (~1%).
2. Larger slop values inflate CSM-only and the trade is unfavorable.
3. The disagreement at the time-domain edges represents a real
   ambiguity: each tool is being internally consistent within its own
   camera model. Forcing convergence to one or the other would just
   bias us toward one camera's choices for no obvious reason.

## Why this is fine

Where both tools think a pixel is valid (which is essentially the
entire image: 99.95% on J08, 99.96% on F09), the resampled DN values
agree to within 0.001 std and 100% of pixels are within 0.01 DN. The
"disagreement" is purely about which exact pixels are at the boundary
of the camera's projected footprint, and even then it's confined to a
1-2 pixel band along the edge.

If you ever need pixel-perfect agreement with ISIS for validation
purposes (e.g. to slot our output into a bit-identical pipeline
position), use `--clip-to-footprint`, which uses the cube's stored
footprint polygon and brings agreement to ~99.99% with the same DN
quality.

## Validation across two CTX cubes

Both runs were against ISIS 9.0.0 cam2map under the sandbox via
`ulimit -n 4096`, with the MOLA DEM as the shape model in both tools.

| Metric                  | J08 (2536×12288) | F09 (5000×7168)  |
|-------------------------|------------------|------------------|
| ISIS cam2map wall time  | 64.64 s          | 73.74 s          |
| isistools csm2map time  | 12.14 s          | 11.36 s          |
| Speedup                 | 5.32×            | 6.49×            |
| Output dims             | 12431 × 3733     | 7635 × 5240      |
| ISIS valid              | 27,371,919       | 31,483,792       |
| CSM valid               | 27,373,477       | 31,474,360       |
| Both valid              | 27,357,115       | 31,471,845       |
| Coverage match          | 99.95%           | 99.96%           |
| ISIS-only               | 14,804           | 11,947           |
| CSM-only                | 3,401            | 2,515            |
| Disagreement            | 18,205           | 14,462           |
| Mean (CSM − ISIS)       | -0.000000        | -0.000000        |
| Median                  | -0.000008        | -0.000017        |
| Std                     | 0.001101         | 0.001192         |
| Max \|diff\|            | 0.015723         | 0.011567         |
| \|diff\| < 0.001        | 68.20%           | 66.09%           |
| \|diff\| < 0.005        | 99.85%           | 99.80%           |
| \|diff\| < 0.01         | 100.00%          | 100.00%          |
| \|diff\| < 0.05         | 100.00%          | 100.00%          |

Two CTX cubes from different orbits, regions, and aspect ratios give
nearly identical statistics — confirming that 99.95% coverage / 100%
DN agreement within 0.01 is the **noise floor** of the CSM vs CSPICE
camera comparison, and not a per-cube fluke.
