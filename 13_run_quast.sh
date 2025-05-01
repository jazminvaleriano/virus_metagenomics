#!/bin/bash
#SBATCH --job-name=quast_evaluation
#SBATCH --output=logs/quast_evaluation_%A.out
#SBATCH --error=logs/quast_evaluation_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate conda environment with BLAST
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate quast_env

# Paths 
ASSEMBLIES_DIR=results/flye_assembly
OUT_DIR=results/flye_assembly/evaluation

mkdir -p "$OUT_DIR"

for assembly in "$ASSEMBLIES_DIR"/barcode*/assembly.fasta; do
    barcode=$(basename "$(dirname "$assembly")")  # Get 'barcodeXX'
    output_dir="$OUT_DIR/$barcode"
    mkdir -p "$output_dir"

    echo "Running MetaQUAST on $barcode..."
    
    metaquast.py "$assembly" \
        -o "$output_dir" \
        --threads 8 \
        --max-ref-number 50 \
        --silent
done