#!/bin/bash
#SBATCH --job-name=kraken2_build
#SBATCH --output=logs/bracken_db%A.out
#SBATCH --error=logs/bracken_db%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate conda environment with bracken
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate bracken_env

KRAKEN_DB="databases/kraken2_std"

# Bracken parameters
KMER_LENGTH=35         # Default for Kraken2
READ_LENGTH=800        # Nanopore median
THREADS=32

# Build Bracken DB
bracken-build -d "$KRAKEN_DB" -t "$THREADS" -k "$KMER_LENGTH" -l "$READ_LENGTH"