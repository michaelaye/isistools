#!/bin/bash
# Install all isistools dependencies via conda-forge.
# Run with the target conda env already activated.
#
# Usage:
#   conda activate <your_env>
#   bash scripts/conda_install_deps.sh [--csm] [--dev]

set -euo pipefail

# Core dependencies
CORE=(
    geopandas
    holoviews
    hvplot
    datashader
    panel
    bokeh
    rioxarray
    xarray
    plio
    pvl
    shapely
    numpy
    pandas
    typer
    matplotlib
    mplcursors
    kalasiris
    pyqt6
    diskcache
)

# CSM deps (csm2map)
CSM=(
    usgscsm
    ale
    scipy
    pyproj
    rich
    rasterio
)

# Dev deps
DEV=(
    pytest
    pytest-cov
    ruff
)

PACKAGES=("${CORE[@]}")

for arg in "$@"; do
    case "$arg" in
        --csm) PACKAGES+=("${CSM[@]}") ;;
        --dev) PACKAGES+=("${DEV[@]}") ;;
        *)     echo "Unknown option: $arg"; echo "Usage: $0 [--csm] [--dev]"; exit 1 ;;
    esac
done

echo "Installing ${#PACKAGES[@]} packages into $(conda info --envs | grep '*' | awk '{print $1}')..."
conda install -c conda-forge "${PACKAGES[@]}"
