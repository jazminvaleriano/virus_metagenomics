#!/bin/bash
#SBATCH --job-name=violin_reads
#SBATCH --output=logs/violin_reads_%A.out
#SBATCH --error=logs/violin_reads_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# Make directory for output
mkdir -p figures 

# Run the Python script
python scripts/lengths_violin_plots.py