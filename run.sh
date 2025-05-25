#!/bin/bash
# filepath: install_linux_deps.sh
# Script to install OpenMP and SFML 3.0.0 for SPH simulation on headless Linux clusters

set -e  # Exit on error
echo "===== Installing Dependencies for SPH Simulation (Headless Linux Cluster) ====="

# Create build and deps directories
WORK_DIR=$(pwd)
DEPS_DIR="${WORK_DIR}/deps"
BUILD_DIR="${WORK_DIR}/build"

mkdir -p "${DEPS_DIR}"
mkdir -p "${BUILD_DIR}"

# Determine if we have sudo access
HAS_SUDO=0
if command -v sudo &> /dev/null; then
    sudo -n true &> /dev/null && HAS_SUDO=1
fi

# Install OpenMP and build essentials
echo "===== Installing OpenMP ====="

if [ $HAS_SUDO -eq 1 ]; then
    sudo apt-get update
    sudo apt-get install -y libomp-dev build-essential cmake
    echo "OpenMP installed via apt"
else
    echo "No sudo access detected - assuming dependencies are available on the cluster"
    echo "If build fails, please ask cluster admin to install: libomp-dev build-essential cmake"
fi

# Install minimal SFML dependencies for headless environment
echo "===== Installing SFML 3.0.0 (Headless Mode) ====="

if [ $HAS_SUDO -eq 1 ]; then
    # Only install dependencies needed for SFML system and network components
    sudo apt-get install -y \
        libudev-dev \
        libfreetype-dev \
        libopenal-dev \
        libflac-dev \
        libvorbis-dev
    echo "Minimal SFML dependencies installed via apt"
else
    echo "No sudo access - cannot install SFML dependencies"
    echo "If build fails, ask cluster admin to install minimal SFML dependencies"
fi

# Download and build SFML 3.0.0
cd "${DEPS_DIR}"
echo "Downloading SFML 3.0.0..."
if [[ ! -d "${DEPS_DIR}/SFML" ]]; then
    git clone --branch 3.0.0 --depth 1 https://github.com/SFML/SFML.git
    cd SFML
else
    cd SFML
    git fetch origin
    git checkout 3.0.0
fi

# Configure and build SFML
echo "Building SFML 3.0.0..."
mkdir -p build && cd build

# Build with appropriate CMake options for headless environment
CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=TRUE"
# Disable graphics components for headless cluster
CMAKE_OPTIONS="$CMAKE_OPTIONS -DSFML_BUILD_WINDOW=FALSE -DSFML_BUILD_GRAPHICS=FALSE"

cmake $CMAKE_OPTIONS -DCMAKE_INSTALL_PREFIX="${DEPS_DIR}/SFML-install" ..

# Determine number of cores for parallel build
NUM_CORES=$(nproc)
echo "Building with $NUM_CORES cores..."
make -j$NUM_CORES
make install

echo "SFML 3.0.0 installed to ${DEPS_DIR}/SFML-install"

# Create environment setup script
cd "${WORK_DIR}"
cat > setup_env.sh << EOF
#!/bin/bash
# Source this file to set up environment variables
# Usage: source setup_env.sh

export SFML_DIR="${DEPS_DIR}/SFML-install"
export LD_LIBRARY_PATH="\${SFML_DIR}/lib:\${LD_LIBRARY_PATH}"
export PATH="\${SFML_DIR}/bin:\${PATH}"

# OpenMP settings optimized for HPC environments
export OMP_NUM_THREADS=\${OMP_NUM_THREADS:-$NUM_CORES}
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OMP_SCHEDULE="dynamic,64"
export OMP_STACKSIZE=16M

echo "Environment set up for SPH simulation"
echo "Using SFML from: \${SFML_DIR}"
echo "OpenMP threads: \${OMP_NUM_THREADS}"
EOF

chmod +x setup_env.sh

# Create build script
cat > build_project.sh << EOF
#!/bin/bash
# Build the SPH project
source ./setup_env.sh

cd "${BUILD_DIR}"
cmake -DCMAKE_BUILD_TYPE=Release \
      -DSFML_DIR="\${SFML_DIR}/lib/cmake/SFML" \
      -DCMAKE_CXX_FLAGS="-O3 -march=native -fopenmp" \
      ..

make -j\${OMP_NUM_THREADS}

echo "Build complete. You can run the application with:"
echo "cd ${BUILD_DIR} && ./sph_simulation"
EOF

chmod +x build_project.sh

# Create a run script for cluster environments
cat > run_cluster.sh << EOF
#!/bin/bash
#PBS -N sph_simulation
#PBS -l nodes=1:ppn=${NUM_CORES}
#PBS -l walltime=01:00:00
# Change the above to match your cluster's job scheduler syntax
# Common alternatives:
# SLURM: #SBATCH --ntasks=1 --cpus-per-task=${NUM_CORES} --time=01:00:00
# LSF: #BSUB -n ${NUM_CORES} -W 60

# Setup environment
cd \${PBS_O_WORKDIR:-\$(pwd)}
source ./setup_env.sh

# Run simulation with performance tracking
cd "${BUILD_DIR}"
echo "Starting SPH simulation with \${OMP_NUM_THREADS} threads at \$(date)"
time ./sph_simulation

echo "Simulation completed at \$(date)"
EOF

chmod +x run_cluster.sh

echo "===== Installation Complete ====="
echo ""
echo "To set up your environment, run:"
echo "  source ./setup_env.sh"
echo ""
echo "To build your project, run:"
echo "  ./build_project.sh"
echo ""
echo "To submit a cluster job (adjust for your scheduler):"
echo "  qsub run_cluster.sh   # PBS"
echo "  sbatch run_cluster.sh # SLURM"
echo "  bsub < run_cluster.sh # LSF"
echo ""
echo "Dependencies installed at: ${DEPS_DIR}"