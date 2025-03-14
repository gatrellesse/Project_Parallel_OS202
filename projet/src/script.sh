#!/bin/bash

CSV_FILE="results.csv"
ENV_VAR_NAME="OMP_NUM_THREADS"   # Replace with your env variable name
EXECUTABLE_PATH="./simulation.exe"  # Replace with your executable path
ENV_VALUES=(1 2 3 4 6 7 8)  # Replace with your test values

# Create CSV header
echo "ENV_VALUE,n_iterations,global_time,avg_time_update,avg_time_affichage" > "$CSV_FILE"

for val in "${ENV_VALUES[@]}"; do
  # Run executable with environment variable and capture full output
  full_output=$(export $ENV_VAR_NAME=$val && $EXECUTABLE_PATH)

  # Parse metrics from output (using awk to extract numbers)
  n_iter=$(echo "$full_output" | awk '/n_iterations/ {print $2}')
  global_t=$(echo "$full_output" | awk '/global_time/ {print $2}')
  avg_time_update=$(echo "$full_output" | awk '/avg_time_update/ {print $2}')
  avg_time_affichage=$(echo "$full_output" | awk '/avg_time_affichage/ {print $2}')

  # Append to CSV
  echo "$val,$n_iter,$global_t,$avg_time_update,$avg_time_affichage" >> "$CSV_FILE"
done

echo "CSV generated: $CSV_FILE"
