#!/bin/bash
set -e

cd /Volumes/Meida/01-CodeLab/sph-simd && 
rm -rf build_headless && 
mkdir build_headless && 
cd build_headless && 
cmake -G Ninja -DHEADLESS_MODE=ON .. && 
ninja && 
./sph_simulation 8000 20 && 
cd ..