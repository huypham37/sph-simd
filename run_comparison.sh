#!/bin/bash

echo "====== SPH SIMD vs Non-SIMD Performance Comparison ======"
echo "Date: $(date)"
echo "System: $(uname -a)"
echo "==========================================================="

# Test duration in seconds
DURATION=60

# Function to run test and extract metrics
run_test() {
    local version="$1"
    local executable="$2"
    local output_file="results_${version}_$(date +%Y%m%d_%H%M%S).log"
    
    echo ""
    echo "Testing $version for ${DURATION} seconds..."
    echo "Executable: $executable"
    
    # Run simulation and capture output
    "$executable" > "$output_file" 2>&1 &
    local sim_pid=$!
    
    # Wait for specified duration
    sleep $DURATION
    
    # Kill the simulation
    kill $sim_pid 2>/dev/null
    wait $sim_pid 2>/dev/null
    
    # Extract metrics from output
    local particles_per_sec=$(grep "Particles/second:" "$output_file" | tail -1 | awk '{print $2}')
    local timesteps_per_sec=$(grep "Timesteps/second:" "$output_file" | tail -1 | awk '{print $2}')
    local avg_time_per_step=$(grep "Avg time per step:" "$output_file" | tail -1 | awk '{print $5}')
    
    echo "Results for $version:"
    echo "  Particles/second: $particles_per_sec"
    echo "  Timesteps/second: $timesteps_per_sec"
    echo "  Avg time per step: $avg_time_per_step"
    
    # Store results in variables for comparison
    if [ "$version" = "SIMD" ]; then
        SIMD_PARTICLES_SEC="$particles_per_sec"
        SIMD_TIMESTEPS_SEC="$timesteps_per_sec"
        SIMD_AVG_TIME="$avg_time_per_step"
    else
        NO_SIMD_PARTICLES_SEC="$particles_per_sec"
        NO_SIMD_TIMESTEPS_SEC="$timesteps_per_sec"
        NO_SIMD_AVG_TIME="$avg_time_per_step"
    fi
    
    echo "Log saved to: $output_file"
}

# Run tests
run_test "SIMD" "./build/sph_simulation"
run_test "NO_SIMD" "./build_no_simd/sph_simulation"

# Calculate improvements
echo ""
echo "==========================================================="
echo "                  PERFORMANCE COMPARISON                   "
echo "==========================================================="

if [ -n "$SIMD_PARTICLES_SEC" ] && [ -n "$NO_SIMD_PARTICLES_SEC" ]; then
    # Calculate percentage improvements
    particles_improvement=$(echo "scale=2; ($SIMD_PARTICLES_SEC - $NO_SIMD_PARTICLES_SEC) / $NO_SIMD_PARTICLES_SEC * 100" | bc -l)
    timesteps_improvement=$(echo "scale=2; ($SIMD_TIMESTEPS_SEC - $NO_SIMD_TIMESTEPS_SEC) / $NO_SIMD_TIMESTEPS_SEC * 100" | bc -l)
    
    echo "Metric                    | NO_SIMD      | SIMD         | Improvement"
    echo "========================= | ============ | ============ | ==========="
    printf "Particles/second          | %-12s | %-12s | +%.1f%%\n" "$NO_SIMD_PARTICLES_SEC" "$SIMD_PARTICLES_SEC" "$particles_improvement"
    printf "Timesteps/second          | %-12s | %-12s | +%.1f%%\n" "$NO_SIMD_TIMESTEPS_SEC" "$SIMD_TIMESTEPS_SEC" "$timesteps_improvement"
    printf "Avg time per step (ms)    | %-12s | %-12s |\n" "$NO_SIMD_AVG_TIME" "$SIMD_AVG_TIME"
    
    echo ""
    echo "SIMD Speedup: ${particles_improvement}% faster particle processing"
    echo "SIMD Speedup: ${timesteps_improvement}% faster simulation stepping"
else
    echo "Error: Could not extract performance metrics from both tests"
fi

echo "==========================================================="
