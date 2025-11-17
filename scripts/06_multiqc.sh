#!/bin/bash

source ../config/config.sh

# Define file paths.
input_dir="$proj_dir_win/data"
output_dir="$proj_dir_win/results/multiqc"
yaml_file="$proj_dir_win/config/multiqc_config.yaml"

# Run MultiQC.
multiqc "$input_dir" \
       -o "$output_dir" \
       -c "$yaml_file"