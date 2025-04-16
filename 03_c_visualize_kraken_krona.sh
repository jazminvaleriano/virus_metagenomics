#!/bin/bash
#SBATCH --job-name=kraken_krona
#SBATCH --output=logs/kraken_krona%A.out
#SBATCH --error=logs/kraken_krona%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate Conda environment (with Krona installed)
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuge_env

# Directory with Kraken output
INPUT_DIR=results/QC/kraken_targeted_reads
OUTPUT_DIR=results/QC/kraken_targeted_reads/krona
mkdir -p "$OUTPUT_DIR"

# Loop over all .kraken files
for FILE in "$INPUT_DIR"/*.kraken2.out; do
    BASENAME=$(basename "$FILE" .kraken2.out)
    OUT_HTML="$OUTPUT_DIR/${BASENAME}.kraken.krona.html"

    echo "Generating Krona plot for $BASENAME..."
    ktImportTaxonomy -o "$OUT_HTML" "$FILE"
done

conda deactivate
