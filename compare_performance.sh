#!/bin/bash

# SPH SIMD vs Non-SIMD Performance Comparison Script

echo "====== SPH SIMD Performance Comparison ======"
echo "Testing SIMD vs Non-SIMD performance"
echo "Date: $(date)"
echo "System: $(uname -a)"
echo "=============================================="

# Create results directory
mkdir -p comparison_results

# Function to build with specific flags
build_version() {
    local version_name="$1"
    local cmake_flags="$2"
    local build_dir="build_${version_name}"
    
    echo "Building $version_name version..."
    echo "CMake flags: $cmake_flags"
    
    # Clean and create build directory
    rm -rf "$build_dir"
    mkdir -p "$build_dir"
    cd "$build_dir"
    
    # Configure and build
    cmake .. $cmake_flags || (echo "CMake failed for $version_name" && cd .. && return 1)
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu) || (echo "Make failed for $version_name" && cd .. && return 1)
    
    cd ..
    echo "$version_name build complete"
    echo ""
}

# Function to run performance test
run_performance_test() {
    local version_name="$1"
    local executable_path="$2"
    local duration="$3"
    local output_file="comparison_results/${version_name}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Testing $version_name performance for ${duration} seconds..."
    
    # Run simulation in background and capture output
    "$executable_path" > "$output_file" 2>&1 &
    local sim_pid=$!
    
    # Wait for specified duration
    sleep $duration
    
    # Stop the simulation
    kill $sim_pid 2>/dev/null
    wait $sim_pid 2>/dev/null
    
    # Extract and display performance metrics
    echo "Results for $version_name:"
    echo "------------------------"
    local metrics=$(tail -20 "$output_file" | grep -A 10 "Performance Metrics" | tail -10)
    echo "$metrics"
    echo "------------------------"
    echo ""
    
    # Extract key metrics for comparison
    local particles_per_sec=$(echo "$metrics" | grep "Particles/second" | awk '{print $2}')
    local timesteps_per_sec=$(echo "$metrics" | grep "Timesteps/second" | awk '{print $2}')
    local avg_time_per_step=$(echo "$metrics" | grep "Avg time per step" | awk '{print $5}')
    
    echo "$version_name,$particles_per_sec,$timesteps_per_sec,$avg_time_per_step" >> comparison_results/summary.csv
}

# Initialize summary CSV
echo "Version,Particles/sec,Timesteps/sec,Avg_time_per_step_ms" > comparison_results/summary.csv

# Build SIMD optimized version (current)
echo "Using existing SIMD-optimized build..."

# Build non-SIMD version
build_version "no_simd" "-DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-O3\ -fno-vectorize"

# Run performance tests
echo "Running performance comparison tests..."
echo ""

# Test SIMD version
run_performance_test "SIMD_Optimized" "./build/sph_simulation" 60

# Test non-SIMD version
run_performance_test "No_SIMD" "./build_no_simd/sph_simulation" 60

# Display comparison summary
echo "====== Performance Comparison Summary ======"
echo ""
cat comparison_results/summary.csv | column -t -s ','
echo ""

# Calculate performance improvement
if [ -f comparison_results/summary.csv ]; then
    simd_timesteps=$(tail -1 comparison_results/summary.csv | cut -d',' -f3)
    no_simd_timesteps=$(head -2 comparison_results/summary.csv | tail -1 | cut -d',' -f3)
    
    if [ ! -z "$simd_timesteps" ] && [ ! -z "$no_simd_timesteps" ]; then
        improvement=$(echo "scale=2; ($simd_timesteps - $no_simd_timesteps) / $no_simd_timesteps * 100" | bc -l 2>/dev/null || echo "calculation not available")
        echo "SIMD Performance Improvement: ${improvement}%"
    fi
fi

echo "=============================================="
