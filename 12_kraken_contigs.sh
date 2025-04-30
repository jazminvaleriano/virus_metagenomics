#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --output=logs/kraken2_%A.out
#SBATCH --error=logs/kraken2_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Kraken2  # or conda activate kraken2_env

DB_DIR="databases/kraken2_custom"

#READS_DIR="filtered_reads" 
READS_DIR="reads_targeted" # Used as QC step

#OUT_DIR="results/kraken2"
OUT_DIR="results/kraken2QC"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq; do #fastq.gz if compressed
    sample=$(basename "$fq" .fastq)

    kraken2 \
        --db "$DB_DIR" \
        --threads 4 \
        --report "$OUT_DIR/${sample}.report.txt" \
        --output "$OUT_DIR/${sample}.kraken2.out" \
        "$fq" \
        #--gzip-compressed 
        
done
