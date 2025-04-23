#!/bin/bash

# Check if both input file and index are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file index"
    exit 1
fi

input_file="$1"
index="$2"
output_file="id_line_uniq_elements.txt"

if [ $index -eq 0 ]; then
    echo ">>> Whole headers" >> "$output_file"
    grep "^@" "$input_file" | sort | uniq >> "$output_file"
else 
    # Find unique elements at the specified position after splitting at colon
    echo ">>> unique ${index}" >> "$output_file"
    grep "^@" "$input_file" | awk -F ":" '{print $'"$index"'}' | sort | uniq >> "$output_file"
fi