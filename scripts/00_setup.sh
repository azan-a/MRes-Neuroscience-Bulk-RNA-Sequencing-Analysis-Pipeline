#!/bin/bash

source ../config/config.sh

# Create required folders on Linux
echo "Setting up project directory on Linux ..."

mkdir -p ~/rna_seq_analysis/raw_fastq
mkdir -p ~/rna_seq_analysis/trimmed_fastq
mkdir -p ~/rna_seq_analysis/fastqc_report

# Create required folders on Windows
echo "Setting up project directory on Windows ..."

# Create all /data/raw folders
mkdir -p "$proj_dir/data/raw/fastq"
mkdir -p "$proj_dir/data/raw/fastqc_reports"
mkdir -p "$proj_dir/data/raw/reference"

# Create all /data/processed folders
mkdir -p "$proj_dir/data/processed/trimmed_fastqc"
mkdir -p "$proj_dir/data/processed/trimmed_fastqc_reports"
mkdir -p "$proj_dir/data/processed/salmon_index"
mkdir -p "$proj_dir/data/processed/salmon_quant"
mkdir -p "$proj_dir/data/processed/rds"
mkdir -p "$proj_dir/data/processed/pathview"

# Create all /results folders
mkdir -p "$proj_dir/results/figures/multifactor"
mkdir -p "$proj_dir/results/figures/singlefactor"
mkdir -p "$proj_dir/results/qc_plots/multifactor"
mkdir -p "$proj_dir/results/qc_plots/singlefactor"
mkdir -p "$proj_dir/results/tables/deg"
mkdir -p "$proj_dir/results/tables/ora"
mkdir -p "$proj_dir/results/tables/gsea"
mkdir -p "$proj_dir/results/sessioninfo"
mkdir -p "$proj_dir/results/multiqc"











