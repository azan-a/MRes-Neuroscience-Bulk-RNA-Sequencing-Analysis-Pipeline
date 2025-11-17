#!/bin/bash

source ../config/config.sh

# Check if only one argument was passed.
if [ $# -ne 1 ]; then
  echo "Invalid argument. Use 'raw' or 'trimmed'"
  exit 1
fi


trim_status="$1"

# Set paths.
if [ "$trim_status" = "raw" ]; then
  echo Analysing raw files
  input_win="$proj_dir_win/data/raw/fastq"
  input_wsl=~/rna_seq_analysis/raw_fastq
  output_win="$proj_dir_win/data/raw/fastqc_reports"
  output_wsl=~/rna_seq_analysis/fastqc_report

elif [ "$trim_status" = "trimmed" ]; then
  echo Analysing trimmed files
  input_win="$proj_dir_win/data/processed/trimmed_fastq"
  input_wsl=~/rna_seq_analysis/trimmed_fastq
  output_win="$proj_dir_win/data/processed/trimmed_fastqc_reports"
  output_wsl=~/rna_seq_analysis/fastqc_report

else
  echo "Invalid argument. Use 'raw' or 'trimmed'"
  exit 1
fi



# Check if input directory in WSL contains files.
if find "$input_wsl" -mindepth 0 -maxdepth 0 -empty | read; then
  echo "Copying fastq files to $input_wsl"
  cp "$input_win"/*.fastq.gz "$input_wsl"/
fi


# Run FastQC
for file in "$input_wsl"/*.fastq.gz; do
  fastqc -o "$output_wsl" "$file"
done


# Move output files to output directory in Windows.
echo "Moving FastQC report to $output_win"
mv "$output_wsl"/*_fastqc.html "$output_win"/
mv "$output_wsl"/*_fastqc.zip "$output_win"/

echo "Move completed."