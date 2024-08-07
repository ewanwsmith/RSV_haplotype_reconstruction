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
echo "n_haps,BIC" > "$output_file"

# Loop through all .out files in the specified directory
for file in "$folder_path"/*.out; do
  # Check if the file exists
  if [ -f "$file" ]; then
    # Read the first line of the file
    first_line=$(head -n 1 "$file")
    # Extract the number after "Final BIC" using awk
    bic_number=$(echo "$first_line" | awk '{for (i=1; i<=NF; i++) if ($i=="BIC") print $(i+1)}')
    # Get the base name of the file (without path)
    file_name=$(basename "$file")
    # Extract the run number from the filename
    n_haps=$(echo "$file_name" | awk -F'_' '{print $2}')
    # Check if the number was found
    if [ -n "$bic_number" ]; then
      # Write the n_haps number and BIC number to the CSV file
      echo "$n_haps,$bic_number" >> "$output_file"
    else
      # Write the n_haps number and "N/A" to the CSV file if no number was found
      echo "$n_haps,N/A" >> "$output_file"
    fi
  fi
done

echo "Extraction complete. Results saved in $output_file."

# Check if bic_plot.R exists in the folder
if [ -f "$folder_path/bic_plot.R" ]; then
  echo "bic_plot.R found. Running the script."
  Rscript "$folder_path/bic_plot.R" "$output_file"
else
  echo "bic_plot.R not found. Skipping the R script execution."
fi