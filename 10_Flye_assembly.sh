#!/bin/bash
#SBATCH --job-name=flye_assembly
#SBATCH --output=logs/flye_assembly_%A.out
#SBATCH --error=logs/flye_assembly_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load Flye 
module load Flye

# Set your input/output
READS="metagenomics_project/filtered_reads"
OUTDIR="results/flye_assembly"

for fq in "$READS"/*.fastq.gz; do
    sample=$(basename "$fq" _filtered.fastq.gz)

    # Run Flye in metagenome mode
    flye --nano-raw "$fq" \
        --out-dir "$OUTDIR" \
        --meta \
        --threads $SLURM_CPUS_PER_TASK