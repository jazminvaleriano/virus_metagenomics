#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH --output=/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/logs/nanoplot_%A.out
#SBATCH --error=/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/logs/nanoplot_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate nanoplot38 # Used a dedicated environment as it needs older python version (3.8)

WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/00_Preprocessing
INPUT_DIR=$WORKDIR/01_length-filtered_reads
OUT_DIR=$WORKDIR/01_QC_length-filt_reads

mkdir -p "$OUT_DIR"

for fq in "$INPUT_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)
    mkdir -p "$OUT_DIR/$sample"
    NanoPlot --fastq "$fq" --loglength --N50 --threads 4 -o "$OUT_DIR/$sample"
done

conda deactivate
