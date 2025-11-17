#!/bin/bash

source ../config/config.sh

# Tell the script where it can find the following folders.
input_win="$proj_dir_win/data/raw/fastq"
input_wsl=~/rna_seq_analysis/raw_fastq
output_win="$proj_dir_win/data/processed/trimmed_fastq"
output_wsl=~/rna_seq_analysis/trimmed_fastq


# Check if input directory in WSL is empty.
if find "$input_wsl" -mindepth 0 -maxdepth 0 -empty | read; then
  echo "Copying fastq files to $input_wsl"
  mv "$input_win"/*.fastq.gz "$input_wsl"/
fi


threads=12


for file in "$input_wsl"/*.gz; do

  # Names the output file to basename + _trimmed.fastq.gz.
  base_name=$(basename "$file" .fastq.gz)
  output_file="$output_wsl/${base_name}_trimmed.fastq.gz"

  # Announce what file is being analysed.
  echo Currently processing $file
  
  # Executes the trimmomatic command.
  trimmomatic SE -threads $threads -phred33 "$file" \
  "$output_file" ILLUMINACLIP:"$adapter_file":2:30:8 \
  2> "${output_file%_trimmed.fastq.gz}_trimmomatic.log"

  echo "Processed $file"
done

echo "All files processed :D"

# Copy output files to output directory in Windows.
echo "Copying trimmed fastq files to $output_win"
cp "$output_wsl"/*.fastq.gz "$output_win"/
mv "$output_wsl"/*.log "$output_win"/

echo "Move completed."
