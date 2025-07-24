#!/bin/bash

## run in interactive session

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# List of accessions to download
accessions=(
  ERR12157771
  ERR12157772
  ERR12157773
)

# Create and enter output directory
mkdir -p 00_reads
cd 00_reads || exit

# Download and convert each accession
for acc in "${accessions[@]}"; do
#   echo "Downloading $acc with prefetch..."
#   prefetch "$acc"

  echo "Converting $acc to FASTQ with fasterq-dump..."
  fasterq-dump "$acc" --split-files

  echo "Done with $acc"
done

echo "All downloads complete."

conda deactivate 