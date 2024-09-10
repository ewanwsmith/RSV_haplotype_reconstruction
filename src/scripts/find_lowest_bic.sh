#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 -in file1.out file2.out ... fileN.out"
    exit 1
}

# Check if the -in flag is provided
if [ "$1" != "-in" ]; then
    usage
fi

shift # Shift to remove the -in argument

# Check if there are files provided
if [ $# -lt 1 ]; then
    usage
fi

# Initialize variables
lowest_bic_file=""
lowest_bic_value=""

# Loop through each provided .out file
for file in "$@"; do
    # Extract the BIC value from the first line of the file
    bic_value=$(awk 'NR==1 {print $3}' "$file")

    # Check if this is the first file or if the current BIC value is lower than the previous lowest
    if [ -z "$lowest_bic_value" ] || (( $(echo "$bic_value < $lowest_bic_value" | bc -l) )); then
        lowest_bic_value="$bic_value"
        lowest_bic_file="$file"
    fi
done

# Output the file with the lowest BIC value
echo "$lowest_bic_file"
