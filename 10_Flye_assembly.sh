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
module load Flye/2.9-GCC-10.3.0

# Set your input/output
READS="filtered_reads"
OUTDIR="results/flye_assembly_hq"

for fq in "$READS"/*_filtered.fastq.gz; do
    sample=$(basename "$fq" _filtered.fastq.gz)
    mkdir -p $OUTDIR/$sample

    # Run Flye in metagenome mode
    flye --nano-hq "$fq" \
        --out-dir $OUTDIR/$sample \
        --meta \
        --threads $SLURM_CPUS_PER_TASK

done