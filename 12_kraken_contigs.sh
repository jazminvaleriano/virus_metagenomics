#!/bin/bash
#SBATCH --job-name=kraken2_contigs
#SBATCH --output=logs/kraken2_contigs_%A.out
#SBATCH --error=logs/kraken2_contigs_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load Kraken2 module 
module load Kraken2

# Path to Kraken2 DB
DB_DIR="databases/kraken2_custom"

# Input directory with contigs (from Flye assemblies)
CONTIGS_DIR="results/flye_assembly"

# Output directory
OUT_DIR="results/kraken2_contigs"
mkdir -p "$OUT_DIR"

# Loop over each barcode assembly file
for fasta in "$CONTIGS_DIR"/barcode*/assembly.fasta; do
    barcode=$(basename "$(dirname "$fasta")")  # extract 'barcodeXX'

    echo "Running Kraken2 on $barcode contigs..."

    kraken2 \
        --db "$DB_DIR" \
        --threads 4 \
        --report "$OUT_DIR/${barcode}.report.txt" \
        --output "$OUT_DIR/${barcode}.kraken2.out" \
        --use-names \
        "$fasta"
done
