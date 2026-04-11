#!/bin/bash
# Benchmark csm2map vs ISIS cam2map on a CTX image at native resolution.
#
# Usage:
#   conda activate isis
#   bash scripts/benchmark_csm2map.sh <input.cub> [resolution_mpp]
#
# Defaults to 6 m/pix if resolution is omitted.

set -euo pipefail

# Sandbox fix — see scripts/isis_sandbox_fix.md
ulimit -n 4096

CUBE="${1:?usage: $0 <input.cub> [resolution_mpp]}"
RESOLUTION="${2:-6.0}"

BASENAME=$(basename "$CUBE" .lev1.cub)
OUTDIR="/tmp/csm2map_bench_${BASENAME}_${RESOLUTION}mpp"
mkdir -p "$OUTDIR"

ISIS_OUTPUT="$OUTDIR/isis_cam2map.cub"
CSM_OUTPUT="$OUTDIR/csm2map.tif"
MAP_FILE="$OUTDIR/equirectangular.map"

CUBE_SIZE=$(du -h "$CUBE" | cut -f1)

echo "============================================================"
echo "csm2map vs cam2map benchmark"
echo "  Cube:       $(basename "$CUBE")  ($CUBE_SIZE)"
echo "  Resolution: ${RESOLUTION} m/pix"
echo "  Output dir: $OUTDIR"
echo "============================================================"

# ---------- Ground range from cube ----------
# camrange prints a PVL document to stdout; parse it with the pvl library
# (robust to whitespace/quote changes in ISIS output) rather than grep/awk.
echo ""
echo "[1/5] camrange..."
CAMRANGE_OUT=$(camrange from="$CUBE" 2>&1) || {
    echo "ERROR: camrange failed. Are you in the ISIS conda env?"
    echo "$CAMRANGE_OUT"
    exit 1
}
echo "$CAMRANGE_OUT" > "$OUTDIR/camrange.pvl"

read -r MINLAT MAXLAT MINLON MAXLON CENTERLAT CENTERLON < <(
    CAMRANGE_TEXT="$CAMRANGE_OUT" python3 << 'PY'
import os
import sys
import pvl

label = pvl.loads(os.environ["CAMRANGE_TEXT"])
# camrange emits a UniversalGroundRange group (planetocentric, +east, 360).
group = None
for key, value in label:
    if hasattr(value, "get") and "MinimumLatitude" in value and key == "UniversalGroundRange":
        group = value
        break
if group is None:
    sys.exit("no UniversalGroundRange group in camrange output")

minlat = float(group["MinimumLatitude"])
maxlat = float(group["MaximumLatitude"])
minlon = float(group["MinimumLongitude"])
maxlon = float(group["MaximumLongitude"])
print(minlat, maxlat, minlon, maxlon, (minlat + maxlat) / 2, (minlon + maxlon) / 2)
PY
)

echo "  Lat: $MINLAT to $MAXLAT"
echo "  Lon: $MINLON to $MAXLON"

cat > "$MAP_FILE" << MAPEOF
Group = Mapping
  ProjectionName     = Equirectangular
  CenterLatitude     = $CENTERLAT
  CenterLongitude    = $CENTERLON
  TargetName         = Mars
  EquatorialRadius   = 3396190.0 <meters>
  PolarRadius        = 3376200.0 <meters>
  LatitudeType       = Planetocentric
  LongitudeDirection = PositiveEast
  LongitudeDomain    = 360
  MinimumLatitude    = $MINLAT
  MaximumLatitude    = $MAXLAT
  MinimumLongitude   = $MINLON
  MaximumLongitude   = $MAXLON
  PixelResolution    = $RESOLUTION <meters/pixel>
End_Group
MAPEOF

# ---------- ISIS cam2map ----------
echo ""
echo "[2/5] ISIS cam2map..."
rm -f "$ISIS_OUTPUT"
ISIS_T0=$(python3 -c "import time; print(time.time())")
cam2map from="$CUBE" to="$ISIS_OUTPUT" map="$MAP_FILE" pixres=map > "$OUTDIR/isis_cam2map.log" 2>&1
ISIS_T1=$(python3 -c "import time; print(time.time())")
ISIS_ELAPSED=$(python3 -c "print(f'{$ISIS_T1 - $ISIS_T0:.2f}')")
echo "  ISIS cam2map: ${ISIS_ELAPSED} s"
echo "  Output size:  $(du -h "$ISIS_OUTPUT" | cut -f1)"

# ---------- csm2map ----------
echo ""
echo "[3/5] csm2map (with --profile)..."
rm -f "$CSM_OUTPUT"
CSM_T0=$(python3 -c "import time; print(time.time())")
isistools csm2map "$CUBE" "$CSM_OUTPUT" --map "$MAP_FILE" --profile 2>&1 | tee "$OUTDIR/csm2map.log"
CSM_T1=$(python3 -c "import time; print(time.time())")
CSM_ELAPSED=$(python3 -c "print(f'{$CSM_T1 - $CSM_T0:.2f}')")
echo "  csm2map total: ${CSM_ELAPSED} s"
echo "  Output size:   $(du -h "$CSM_OUTPUT" | cut -f1)"

# ---------- Compare ----------
echo ""
echo "[4/5] Comparing outputs..."
isistools csm2map-compare "$ISIS_OUTPUT" "$CSM_OUTPUT" 2>&1 | tee "$OUTDIR/compare.log"

# ---------- Summary ----------
echo ""
echo "[5/5] Summary:"
echo "============================================================"
printf "  %-20s %10s\n" "cam2map (ISIS):"  "${ISIS_ELAPSED} s"
printf "  %-20s %10s\n" "csm2map (ours):"  "${CSM_ELAPSED} s"
SPEEDUP=$(python3 -c "print(f'{$ISIS_ELAPSED / $CSM_ELAPSED:.2f}x')")
printf "  %-20s %10s\n" "speedup:" "$SPEEDUP"
echo "============================================================"
echo "Logs + outputs in $OUTDIR"
