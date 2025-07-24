#!/bin/bash
#SBATCH --job-name=nanoplot
#SBATCH --output=logs/01_nanoplot.out
#SBATCH --error=logs/01_nanoplot.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate nanoplot38

RAW_DIR="/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/full_reads_sim_om/reads"
OUT_DIR="/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/01_QC_results"

mkdir -p "$OUT_DIR"

for fq in $RAW_DIR/*om.fastq ; do
    sample=$(basename "$fq")
    mkdir -p "$OUT_DIR/$sample"
    ~/.conda/envs/nanoplot38/bin/NanoPlot --fastq "$fq" --loglength --N50 --threads 4 -o "$OUT_DIR/$sample"
done

conda deactivate
