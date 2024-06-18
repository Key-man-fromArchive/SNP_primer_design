#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 -i input_dir -o output_dir"
  exit 1
}

# Parse command line arguments
while getopts "i:o:" opt; do
  case $opt in
    i) input_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    *) usage ;;
  esac
done

# Check if input and output directories are provided
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
  usage
fi

# Define the merged output file name
merged_output_file="$output_dir/merged_output.txt"
merged_fasta_file="$output_dir/Merged.fasta"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all files in the input directory
for input_file in "$input_dir"/*; do
  # Get the base name of the file (without path and extension)
  base_name=$(basename "$input_file" .fasta)
  
  # Define the output file name
  output_file="$output_dir/${base_name}_output.txt"
  
  # Run the Python script and redirect the output to the output file
  python ./ASQprimer-line.py "$input_file" > "$output_file"
done

# Merge all output files into one
cat "$output_dir"/*_output.txt > "$merged_output_file"

# Extract FASTA sequences and save to Merged.fasta
awk '/^>/ {if (seq) print seq; print; seq=""} /^[ACGTNacgtn]+$/ {seq=seq$0} END {if (seq) print seq}' "$merged_output_file" > "$merged_fasta_file"

echo "All outputs have been merged into $merged_output_file"
echo "FASTA sequences have been extracted to $merged_fasta_file"
