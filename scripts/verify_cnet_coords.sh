#!/bin/bash
# Verify cnet XYZ-to-lonlat conversion matches campt output,
# and test footprintinit with MAXEMISSION=89.5 for tighter footprint boundaries.
#
# Run in the isis conda environment:
#   conda activate isis && bash scripts/verify_cnet_coords.sh

CUBE="/Users/maye/planetarypy_data/missions/mro/ctx/mrox_3111/K05_055343_1840_XI_04N287W/K05_055343_1840_XI_04N287W.lev1.cub"
SAMPLE=1407.0
LINE=371.9

echo "=== 1. Verify XYZ-to-lonlat conversion with campt ==="
echo "Expected from XYZ conversion: lon=72.718772, lat_centric=3.485787"
echo ""
campt from="$CUBE" sample="$SAMPLE" line="$LINE" 2>&1 | \
  grep -E 'PlanetocentricLatitude|PlanetographicLatitude|PositiveEast360|BodyFixedCoordinate'

echo ""
echo "=== 2. Re-run footprintinit with MAXEMISSION=89.5 ==="
echo "Current footprint lat range: 3.319 to 4.833"
echo ""
footprintinit from="$CUBE" maxemission=89.5
echo "Exit code: $?"

echo ""
echo "=== 3. Check new footprint bounds ==="
python3 -c "
from isistools.io.footprints import read_footprint
from pathlib import Path
geom = read_footprint(Path('$CUBE'))
bounds = geom.bounds
print(f'New footprint bounds: lon=[{bounds[0]:.6f}, {bounds[2]:.6f}], lat=[{bounds[1]:.6f}, {bounds[3]:.6f}]')
"
