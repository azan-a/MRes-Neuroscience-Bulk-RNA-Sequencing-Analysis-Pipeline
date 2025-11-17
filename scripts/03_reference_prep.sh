#!/bin/bash

source ../config/config.sh

# Input files.
genome_fasta="$proj_dir_win/data/raw/reference/\
GCF_036323735.1_GRCr8_genomic.fna.gz"
transcriptome_fasta="$proj_dir_win/data/raw/reference/\
GCF_036323735.1_GRCr8_rna.fna.gz"

# Output files.
combined_fasta="$proj_dir_win/data/processed/salmon_index/\
transcriptome_plus_genome.fa.gz"
decoys_file="$proj_dir_win/data/processed/salmon_index/decoys.txt"


# Combine RNA and genome FASTA files.
cat "$rna_fasta" "$genome_fasta" > "$combined_fasta"

# Extract genomic headers as decoys.
grep "^>" <(gunzip -c "$genome_fasta") | cut -d " " -f 1 > "$decoys_file"

# Remove leading ">" from decoy names.
sed -i.bak -e 's/>//g' "$decoys_file"

echo Processing done.