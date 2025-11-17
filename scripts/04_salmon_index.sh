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

#Tell the script where it can find the following files.
transcriptome_plus_genome="./transcriptome_plus_genome.fa.gz"
decoys="./decoys.txt"
output_index="./rattus_norvegicus_index/"
salmon_path=

#Create an index for salmon.
$salmon_path index -t $transcriptome_plus_genome -d $decoys -p 2 \
  -i $output_index
