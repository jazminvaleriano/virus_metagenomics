#!/bin/bash
#SBATCH --job-name=asm_evaluation
#SBATCH --output=logs/asm_evaluation_%A.out
#SBATCH --error=logs/asm_evaluation_%A.err
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

mkdir -p figures 

# Run the Python scripts
#python scripts/comparison_tools_GR.py
python scripts/genomes_recovery.py
