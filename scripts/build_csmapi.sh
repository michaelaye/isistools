#!/bin/bash
# Build csmapi Python bindings from source (arm64).
# The conda-forge csmapi package only has x86_64 builds.
#
# Prerequisites: activate a py312 env with the csm C++ library installed:
#   conda activate py312
#   conda install -c conda-forge cmake swig csm
#
# Usage:
#   bash scripts/build_csmapi.sh

set -euo pipefail

BUILDDIR="/tmp/swigcsm"

# Clone if not already present
if [ ! -d "$BUILDDIR" ]; then
    echo "Cloning swigcsm..."
    git clone https://github.com/DOI-USGS/swigcsm.git "$BUILDDIR"
fi

cd "$BUILDDIR"

# Clean previous build
rm -rf build
mkdir build
cd build

echo "Configuring..."
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DPython_EXECUTABLE="$(which python)"

echo "Building..."
cmake --build .

echo "Installing Python package..."
cd python
pip install .

echo ""
echo "Verifying..."
python -c "import csmapi; print('csmapi OK:', csmapi.__file__)"
