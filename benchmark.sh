#!/bin/bash

# SPH SIMD Performance Benchmark Script
# This script runs the simulation multiple times to gather performance statistics

echo "====== SPH SIMD Performance Benchmark ======"
echo "Testing SIMD optimizations on Apple Silicon (M1/M2/M3)"
echo "Date: $(date)"
echo "System: $(uname -a)"
echo "Compiler: $(clang++ --version | head -1)"
echo "=============================================="

# Create results directory
mkdir -p benchmark_results

# Function to run benchmark
run_benchmark() {
    local test_name="$1"
    local duration="$2"
    local output_file="benchmark_results/${test_name}_$(date +%Y%m%d_%H%M%S).log"
    
    echo "Running benchmark: $test_name for ${duration} seconds"
    echo "Output file: $output_file"
    
    # Run simulation in background and capture output
    ./build/sph_simulation > "$output_file" 2>&1 &
    local sim_pid=$!
    
    # Wait for specified duration
    sleep $duration
    
    # Stop the simulation
    kill $sim_pid 2>/dev/null
    wait $sim_pid 2>/dev/null
    
    # Extract performance metrics from the log
    echo "Results for $test_name:"
    echo "------------------------"
    
    # Get the last performance metrics block
    tail -20 "$output_file" | grep -A 10 "Performance Metrics" | tail -10
    echo "------------------------"
    echo ""
}

# Ensure simulation is built
if [ ! -f "./build/sph_simulation" ]; then
    echo "Building simulation..."
    ./build.sh
fi

# Run benchmarks
echo "Starting benchmarks..."
echo ""

# Test 1: Short run for quick metrics
run_benchmark "simd_optimized_quick" 30

# Test 2: Medium run for stable metrics
run_benchmark "simd_optimized_medium" 60

# Test 3: Extended run for thorough testing
run_benchmark "simd_optimized_extended" 120

echo "====== Benchmark Complete ======"
echo "Results saved in benchmark_results/ directory"
echo ""

# Summary of latest results
echo "====== Latest Performance Summary ======"
latest_file=$(ls -t benchmark_results/*.log | head -1)
if [ -f "$latest_file" ]; then
    echo "From: $latest_file"
    tail -15 "$latest_file" | grep -A 10 "Performance Metrics" | tail -10
fi
echo "========================================"
