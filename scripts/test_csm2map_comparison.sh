#!/bin/bash
# Compare ISIS cam2map vs isistools csm2map on a CTX image.
#
# This script runs in two phases:
#   Phase 1 (ISIS):     cam2map + camrange (needs ISIS conda env)
#   Phase 2 (isistools): csm2map + comparison (needs CSM deps)
#
# The ISIS env needs the CSM stack alongside ISIS itself:
#
#   conda activate isis
#   conda install -c conda-forge usgscsm ale rasterio rich
#   pip install -e /Users/maye/Dropbox/src/isistools
#
# Then run:
#   conda activate isis
#   bash scripts/test_csm2map_comparison.sh

set -euo pipefail

# ---------- Sandbox workaround ----------
# Claude Code's sandbox on macOS sets RLIMIT_NOFILE to INT64_MAX, which breaks
# ISIS 9.0.0's CubeManager cache sizing (`p_maxOpenFiles = rlim_cur * 0.60`).
# The overflowed cache limit causes a Qt QList corruption that manifests as
# SIGSEGV in _platform_memmove$VARIANT$Rosetta during DemShape construction.
# Setting a sane file descriptor limit fixes every camera-loading ISIS tool
# (cam2map, campt, camrange, caminfo, etc.) and produces output that is
# bit-identical to the same command run from a normal terminal.
ulimit -n 4096

# ---------- Configuration ----------
CUBE="/Users/maye/planetarypy_data/missions/mro/ctx/mrox_2774/J08_048038_1842_XN_04N287W/J08_048038_1842_XN_04N287W.lev1.cub"
OUTDIR="/tmp/csm2map_comparison"
RESOLUTION=6.0  # meters/pixel

ISIS_OUTPUT="$OUTDIR/J08_isis_cam2map.cub"
CSM_OUTPUT="$OUTDIR/J08_csm2map.tif"
MAP_FILE="$OUTDIR/equirectangular.map"

mkdir -p "$OUTDIR"

echo "============================================================"
echo "csm2map vs ISIS cam2map comparison"
echo "  Cube: $(basename $CUBE)"
echo "  Resolution: ${RESOLUTION} m/pix"
echo "  Output dir: $OUTDIR"
echo "============================================================"

# ---------- Step 1: Get ground range from the cube ----------
echo ""
echo "[Step 1] Extracting ground range with camrange..."
CAMRANGE_OUT=$(camrange from="$CUBE" 2>&1) || {
    echo "ERROR: camrange failed. Are you in the ISIS conda env?"
    exit 1
}

# Parse lat/lon range from camrange output
MINLAT=$(echo "$CAMRANGE_OUT" | grep "MinimumLatitude" | head -1 | awk '{print $3}')
MAXLAT=$(echo "$CAMRANGE_OUT" | grep "MaximumLatitude" | head -1 | awk '{print $3}')
MINLON=$(echo "$CAMRANGE_OUT" | grep "MinimumLongitude" | head -1 | awk '{print $3}')
MAXLON=$(echo "$CAMRANGE_OUT" | grep "MaximumLongitude" | head -1 | awk '{print $3}')

echo "  Lat range: $MINLAT to $MAXLAT"
echo "  Lon range: $MINLON to $MAXLON"

# ---------- Step 2: Create a MAP file with fixed bounds ----------
# Using a MAP file with explicit bounds ensures both tools use identical grids
echo ""
echo "[Step 2] Creating MAP file with explicit bounds..."

# Compute center lat/lon
CENTERLAT=$(python3 -c "print(($MINLAT + $MAXLAT) / 2)")
CENTERLON=$(python3 -c "print(($MINLON + $MAXLON) / 2)")

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

echo "  MAP file: $MAP_FILE"
cat "$MAP_FILE"

# ---------- Step 3: Run ISIS cam2map ----------
echo ""
echo "[Step 3] Running ISIS cam2map..."
rm -f "$ISIS_OUTPUT"
time cam2map from="$CUBE" to="$ISIS_OUTPUT" map="$MAP_FILE" pixres=map
echo "  Output: $ISIS_OUTPUT ($(du -h "$ISIS_OUTPUT" | cut -f1))"

# Print the output mapping group for verification
echo "  Output mapping:"
catlab from="$ISIS_OUTPUT" | grep -A 20 "Group = Mapping" | head -25

# ---------- Step 4: Run isistools csm2map ----------
echo ""
echo "[Step 4] Running isistools csm2map..."
rm -f "$CSM_OUTPUT"
time isistools csm2map "$CUBE" "$CSM_OUTPUT" --map "$MAP_FILE"
echo "  Output: $CSM_OUTPUT ($(du -h "$CSM_OUTPUT" | cut -f1))"

# ---------- Step 5: Compare ----------
echo ""
echo "[Step 5] Comparing outputs..."
isistools csm2map-compare "$ISIS_OUTPUT" "$CSM_OUTPUT"

echo ""
echo "============================================================"
echo "Done. Files in $OUTDIR:"
ls -lh "$OUTDIR"
echo "============================================================"
