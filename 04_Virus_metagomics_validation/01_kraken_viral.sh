#!/bin/bash
#SBATCH --job-name=kraken2_run
#SBATCH --output=logs/kraken2_run%A.out
#SBATCH --error=logs/kraken2_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load Kraken2 module
module load Kraken2/2.1.2-gompi-2021a

# -------------------- CONFIG --------------------

READS_DIR=00_simulated_reads
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kraken2_RefSeqViral
OUT_DIR=1st_classification_results/kraken2
CLASS_READS_DIR=02_virus_enriched_reads/kraken2
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

mkdir -p "$OUT_DIR"
mkdir -p $CLASS_READS_DIR

for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)

    echo "Running Kraken2 on sample: $sample"

    kraken2 \
        --db "$DB_DIR" \
        --threads $THREADS \
        --report "$OUT_DIR/${sample}.report.txt" \
        --output "$OUT_DIR/${sample}.kraken2.out" \
        --unclassified-out "$OUT_DIR/${sample}_unclassified.fastq" \
        --classified-out "$CLASS_READS_DIR/${sample}_classified.fastq" \
        "$fq"

done

echo "Kraken2 run complete."