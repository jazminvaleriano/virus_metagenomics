#!/bin/bash
#SBATCH --job-name=kaiju_krona
#SBATCH --output=logs/kaiju_krona%A.out
#SBATCH --error=logs/kaiju_krona%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate Conda environment (where I have installed Krona)
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate kaiju_env

# Paths to Kaiju taxonomy files
NODES=databases/kaiju/nodes.dmp
NAMES=databases/kaiju/names.dmp
INPUT_DIR=results/kaiju

# Process each .kaiju.out file
for FILE in "$INPUT_DIR"/*.kaiju.out; do
    BASENAME=$(basename "$FILE" .kaiju.out)
    KRONA_TXT="$INPUT_DIR/${BASENAME}.kaiju.krona"
    KRONA_HTML="$INPUT_DIR/${BASENAME}.kaiju.krona.html"

    echo "Processing $BASENAME..."

    kaiju2krona -t "$NODES" -n "$NAMES" -i "$FILE" -o "$KRONA_TXT" -u
    ktImportText -o "$KRONA_HTML" "$KRONA_TXT"
done

conda deactivate
