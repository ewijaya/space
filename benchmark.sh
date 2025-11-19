#!/bin/bash
#
# SPACE Performance Benchmark Script
# Compares performance with different thread counts
#

TESTFILE="files/test.fasta"
SPECIES="SC"

echo "===================================================================="
echo "SPACE Performance Benchmark"
echo "===================================================================="
echo ""
echo "Test file: $TESTFILE"
echo "Species: $SPECIES"
echo ""

# Check if test file exists
if [ ! -f "$TESTFILE" ]; then
    echo "Error: Test file $TESTFILE not found!"
    exit 1
fi

# Get number of CPU cores
NCORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
echo "Available CPU cores: $NCORES"
echo ""

# Make sure binaries are compiled
if [ ! -f "generate_v.exe" ] || [ ! -f "mining_v.exe" ] || [ ! -f "scoring_v.exe" ]; then
    echo "Compiling binaries..."
    ./compileall
    echo ""
fi

# Function to run benchmark
run_benchmark() {
    local threads=$1
    local label=$2

    echo "--------------------------------------------------------------------"
    echo "Running with $label"
    echo "--------------------------------------------------------------------"

    export OMP_NUM_THREADS=$threads

    # Run 3 times and take the average
    local total_time=0
    local runs=3

    for i in $(seq 1 $runs); do
        echo -n "  Run $i/$runs... "

        # Capture timing output
        local output=$(python3 run_all.py $TESTFILE $SPECIES 2>&1)

        # Extract time from each program
        local gen_time=$(echo "$output" | grep -A1 "generate_v.exe" | grep "Total time:" | awk '{print $3}')
        local min_time=$(echo "$output" | grep -A1 "mining_v.exe" | grep "Total time:" | awk '{print $3}')
        local sco_time=$(echo "$output" | grep -A1 "scoring_v.exe" | grep "Total time:" | awk '{print $3}')

        # Calculate total (rough estimate)
        local run_time=$(echo "$gen_time + $min_time + $sco_time" | bc -l 2>/dev/null || echo "N/A")

        if [ "$run_time" != "N/A" ]; then
            echo "${run_time}s (gen: ${gen_time}s, min: ${min_time}s, score: ${sco_time}s)"
            total_time=$(echo "$total_time + $run_time" | bc -l)
        else
            echo "timing data not available"
        fi
    done

    if [ "$total_time" != "0" ]; then
        local avg_time=$(echo "scale=2; $total_time / $runs" | bc -l)
        echo ""
        echo "  Average time: ${avg_time}s"
        echo ""
    fi
}

# Benchmark with different thread counts
echo "===================================================================="
echo "BENCHMARKING WITH DIFFERENT THREAD COUNTS"
echo "===================================================================="
echo ""

# Single thread baseline
run_benchmark 1 "1 thread (baseline)"

# Half cores
if [ $NCORES -gt 2 ]; then
    HALF_CORES=$((NCORES / 2))
    run_benchmark $HALF_CORES "$HALF_CORES threads"
fi

# All cores
if [ $NCORES -gt 1 ]; then
    run_benchmark $NCORES "$NCORES threads (all cores)"
fi

echo "===================================================================="
echo "Benchmark complete!"
echo ""
echo "Note: Speedup is most noticeable in the scoring phase (scoring_v.exe)"
echo "The generate and mining phases are not yet parallelized."
echo "===================================================================="
