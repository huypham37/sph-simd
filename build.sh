#!/bin/bash
set -e

echo "=== SPH Fluid Simulation Build Script ==="

# Make sure build directory exists
mkdir -p build
cd build

# Configure with CMake using Ninja generator
echo "Configuring with CMake and Ninja..."
cmake -G Ninja ..

# Build the project
echo "Building SPH simulation..."
ninja

echo "=== Build complete! ==="
echo "Run the simulation with: ./sph_simulation"

# Make the executable executable just in case
chmod +x ./sph_simulation 2>/dev/null || true