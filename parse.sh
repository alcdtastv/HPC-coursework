#!/bin/bash

# Extract number of threads and real time from input file
input_file="input.txt"
output_file="output.csv"
echo "Threads,Time (s)" > "$output_file"
while read -r line; do
  if [[ "$line" == "Number of threads:"* ]]; then
    threads=${line##*: }
  elif [[ "$line" == "real"* ]]; then
    time=$(echo "$line" | awk -F'\t' '{print $2}' | sed 's/m/:/' | awk -F':' '{print ($1 * 60) + $2}')
    echo "$threads,$time" >> "$output_file"
  fi
done < "$input_file"
