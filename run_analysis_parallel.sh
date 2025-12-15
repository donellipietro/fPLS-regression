#!/bin/bash

# $1: name R script
# $2: name test
# $3: number of calls


start=$(date +%s.%N)

# Define the number of high-performance cores
NUM_CORES=12

# Define the total number of CPU cores
TOTAL_CORES=$(sysctl -n hw.physicalcpu)

# Calculate the number of cores to be used for parallel processing
PARALLEL_CORES=$((TOTAL_CORES - 4))

# Ensure that PARALLEL_CORES does not exceed the available high-performance cores
if [ "$PARALLEL_CORES" -gt "$NUM_CORES" ]; then
  PARALLEL_CORES=$NUM_CORES
fi


###############################################################################


# Run tasks in parallel using GNU Parallel
seq 1 $3 | parallel -j "$PARALLEL_CORES" RScript analysis/$1 $2 {}


###############################################################################


# End timer
end=$(date +%s.%N)
execution_time=$(echo "$end - $start" | bc)

# Print total execution time
echo "Total execution time: $execution_time seconds"
