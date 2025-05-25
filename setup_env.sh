#!/bin/bash
# Source this file to set up environment variables
# Usage: source setup_env.sh

export SFML_DIR="/home/coder/project/sph-simd/deps/SFML-install"
export LD_LIBRARY_PATH="${SFML_DIR}/lib:${LD_LIBRARY_PATH}"
export PATH="${SFML_DIR}/bin:${PATH}"

# OpenMP settings optimized for HPC environments
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-128}
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic,64"
export OMP_STACKSIZE=16M

echo "Environment set up for SPH simulation"
echo "Using SFML from: ${SFML_DIR}"
echo "OpenMP threads: ${OMP_NUM_THREADS}"
