#!/bin/bash

# Function to display usage information
usage() {
    echo "Usage: $0 -i input.fasta -o output.fasta -s start -e end"
    exit 1
}

# Parse command line arguments
while getopts ":i:o:s:e:" opt; do
    case $opt in
        i) input_file="$OPTARG" ;;
        o) output_file="$OPTARG" ;;
        s) start="$OPTARG" ;;
        e) end="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if all arguments are provided
if [[ -z "$input_file" || -z "$output_file" || -z "$start" || -z "$end" ]]; then
    usage
fi

# Ensure start and end are integers and start is less than end
if ! [[ "$start" =~ ^[0-9]+$ ]] || ! [[ "$end" =~ ^[0-9]+$ ]] || [[ "$start" -ge "$end" ]]; then
    echo "Error: Start and end must be integers, and start must be less than end."
    exit 1
fi

# Adjust for zero-based indexing (optional, depending on your requirements)
start=$((start - 1))
end=$((end - 1))

# Process the input FASTA file
{
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            # Header line, print it to the output file
            echo "$line"
        else
            # Sequence line, trim and print to the output file
            echo "${line:$start:$(($end - $start + 1))}"
        fi
    done
} < "$input_file" > "$output_file"

echo "Processing complete. Trimmed sequences saved to $output_file."