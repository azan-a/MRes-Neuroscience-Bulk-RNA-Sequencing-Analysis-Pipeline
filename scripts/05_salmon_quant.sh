#!/bin/bash

#Run with current environment (-V) and in the current directory (-cwd).
#$ -V -cwd

#Request some time- min 15 mins - max 48 hours.
#$ -l h_rt=12:00:00

#Request some memory per core.
#$ -pe smp 2
#$ -l h_vmem=32G

#Get email at start and end of the job.
#$ -m be


# Specify file paths.
input_dir="./trimmed_fastq_files"
output_dir="./salmon_output/"
index_dir="./rattus_norvegicus_index/"
salmon_path=


for file in "$input_dir"/*.gz; do
    # Extract filename without extension for output directory.
    filename=$(basename "$file" .fastq.gz)
    file_output_dir="$output_dir/$filename"
    
    # Create a unique output directory for each input file.
    mkdir -p "$file_output_dir"
    
    # Run Salmon Quantification.
    $salmon_path quant -i "$index_dir" -l A -r "$file" -p 12 \
      -o "$file_output_dir" --seqBias
done
