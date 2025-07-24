#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --output=logs/kraken2_%A.out
#SBATCH --error=logs/kraken2_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Kraken2/2.1.2-gompi-2021a

DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kraken2_PFP
READS_DIR=00_simulated_reads
OUT_DIR=results/01_kraken

mkdir -p "$OUT_DIR"

# Create CSV header
echo "classifier,sample,start_time,end_time,elapsed_seconds" > "runtime_benchmark.csv"


for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)
    
    start_time=$(date +%s)

    echo "[$(date)] Starting classification for $sample"

    kraken2 \
        --db "$DB_DIR" \
        --threads 8 \
        --report "$OUT_DIR/${sample}.report.txt" \
        --output "$OUT_DIR/${sample}.kraken2.out" \
        "$fq"
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    echo "[$(date)] Finished $sample in $elapsed seconds"

    # Log to CSV
    echo "kraken,$sample,$start_time,$end_time,$elapsed" >> "runtime_benchmark.csv"
done



