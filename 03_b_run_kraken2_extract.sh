#!/bin/bash
#SBATCH --job-name=kraken2_unclassified
#SBATCH --output=logs/kraken2_unclassified_%A.out
#SBATCH --error=logs/kraken2_unclassified_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Kraken2/2.1.2-gompi-2021a

DB_DIR="databases/kraken2_bracken"
READS_DIR="filtered_reads/second_basecall"
OUT_DIR="results/kraken_vs_standard"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" _filtered.fastq.gz)

    kraken2 \
        --db "$DB_DIR" \
        --threads 8 \
        --gzip-compressed \
        --report "$OUT_DIR/${sample}.report.txt" \
        --output "$OUT_DIR/${sample}.kraken2.out" \
        --unclassified-out "$OUT_DIR/${sample}_unclassified.fastq" \
        "$fq"

done



