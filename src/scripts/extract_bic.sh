#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <path-to-folder>"
  exit 1
fi

# Get the path to the folder from the argument
folder_path=$1

# Output CSV file
output_file="$folder_path/extracted_bic.csv"

# Write the header to the CSV file
echo "Run,BIC" > "$output_file"

# Loop through all .out files in the specified directory
for file in "$folder_path"/*.out; do
  # Check if the file exists
  if [ -f "$file" ]; then
    # Read the first line of the file
    first_line=$(head -n 1 "$file")
    # Extract the number after "Final BIC"
    bic_number=$(echo "$first_line" | grep -oP 'Final BIC \K\d+')
    # Get the base name of the file (without path)
    file_name=$(basename "$file")
    # Check if the number was found
    if [ -n "$bic_number" ]; then
      # Write the file name and BIC number to the CSV file
      echo "$file_name,$bic_number" >> "$output_file"
    else
      # Write the file name and "N/A" to the CSV file if no number was found
      echo "$file_name,N/A" >> "$output_file"
    fi
  fi
done

echo "Extraction complete. Results saved in $output_file."