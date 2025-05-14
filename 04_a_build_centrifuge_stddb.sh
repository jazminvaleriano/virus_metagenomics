#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH --output=logs/centrifuge_%A.out
#SBATCH --error=logs/centrifuge_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuge_env

DB="databases/centrifuge/p_compressed+h+v"
READS_DIR="filtered_reads"
OUT_DIR="results/centrifuge"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" _filtered.fastq.gz)

    centrifuge \
        -x "$DB" \
        -U "$fq" \
        -S "$OUT_DIR/${sample}.centrifuge.tsv" \
        --report-file "$OUT_DIR/${sample}.centrifuge.report" \
        --threads 4

done

conda deactivate