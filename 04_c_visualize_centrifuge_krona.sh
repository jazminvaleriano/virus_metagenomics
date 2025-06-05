#!/bin/bash
#SBATCH --job-name=centrifuge_krona
#SBATCH --output=logs/centrifuge_krona_%A.out
#SBATCH --error=logs/centrifuge_krona_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuge_env

IN_DIR="results/centrifuge"
DB="databases/centrifuge"

# Set Krona taxonomy path and show debug info
export KRONA_TAXONOMY="$CONDA_PREFIX/opt/krona/taxonomy"
echo "Using taxonomy path: $KRONA_TAXONOMY"
ls "$KRONA_TAXONOMY"

for tsv in "$IN_DIR"/*.centrifuge.tsv; do
    sample=$(basename "$tsv" .centrifuge.tsv)

    # Generate kreport format 
    centrifuge-kreport -x "$DB/p_compressed+h+v" "$tsv" > "$IN_DIR/${sample}.kreport.txt"

    # Generate Krona HTML with explicit taxonomy path
    ktImportTaxonomy "$IN_DIR/${sample}.kreport.txt" \
        -o "$IN_DIR/${sample}.krona.html" \
        -tax "$KRONA_TAXONOMY"
done

conda deactivate