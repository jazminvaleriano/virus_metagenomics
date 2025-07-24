#!/bin/bash
#SBATCH --job-name=kraken2_run
#SBATCH --error=logs/kraken2_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G

# Load Kraken2 module
module load Kraken2/2.1.2-gompi-2021a

# -------------------- CONFIG --------------------

READS_DIR=00_reads
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kraken2_RefSeqViral
OUT_DIR=01_classification_results/kraken
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)

    echo "Running Kraken2 on sample: $sample"

    kraken2 \
        --db "$DB_DIR" \
        --threads $THREADS \
        --report "$OUT_DIR/${sample}.report.txt" \
        --output "$OUT_DIR/${sample}.kraken2.out" \
        --unclassified-out "$OUT_DIR/${sample}_unclassified.fastq" \
        "$fq"

done

echo "Kraken2 run complete."