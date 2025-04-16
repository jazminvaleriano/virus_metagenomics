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
conda activate nanoplot38

RAW_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/raw_data"
OUT_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/qc"

mkdir -p "$OUT_DIR"

for fq in "$RAW_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    mkdir -p "$OUT_DIR/$sample"
    ~/.conda/envs/nanoplot38/bin/NanoPlot --fastq "$fq" --loglength --N50 --threads 4 -o "$OUT_DIR/$sample"
done

conda deactivate
