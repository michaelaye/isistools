#!/bin/bash
# Re-run footprintinit on both cubes after jigsaw to update footprint polygons.
# Run in isis conda env:
#   conda activate isis && bash scripts/rerun_footprintinit.sh

CUBE1="/Users/maye/planetarypy_data/missions/mro/ctx/mrox_3111/K05_055343_1840_XI_04N287W/K05_055343_1840_XI_04N287W.lev1.cub"
CUBE2="/Users/maye/planetarypy_data/missions/mro/ctx/mrox_4882/V05_082075_1840_XN_04N287W/V05_082075_1840_XN_04N287W.lev1.cub"

echo "Re-running footprintinit on both cubes (post-jigsaw)..."

for CUBE in "$CUBE1" "$CUBE2"; do
    echo "  $(basename $CUBE)..."
    footprintinit from="$CUBE" maxemission=89.5
    echo "    exit: $?"
done

echo ""
echo "Checking new footprint bounds..."
python3 -c "
from isistools.io.footprints import load_footprints
from isistools.io.cache import get_cache
from pathlib import Path

# Clear footprint cache so we get fresh polygons
cache = get_cache()
for key in list(cache):
    if key.startswith('footprint:'):
        del cache[key]
print('Cleared footprint cache')

cubelist = Path('/Volumes/planet/Mars/CTX/special/test_isis_gap_pipeline/image_list_lev1.lis')
gdf = load_footprints(cubelist)
print(f'Combined footprint bounds: {gdf.total_bounds}')

from isistools.io.controlnet import load_cnet
from isistools.plotting.cnet_overlay import cnet_to_geodataframe
cnet_df = load_cnet('/Volumes/planet/Mars/CTX/special/test_isis_gap_pipeline/jigged_network.net')
cnet_gdf = cnet_to_geodataframe(cnet_df)
print(f'Cnet point bounds: {cnet_gdf.total_bounds}')

from isistools.plotting.footprint_mpl import footprint_png
out = footprint_png(gdf, '/tmp/footprints_cnet_updated.png', cnet_df=cnet_df, title='Post-jigsaw footprints + cnet')
print(f'Saved: {out}')
"
