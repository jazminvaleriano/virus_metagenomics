#!/bin/bash
#SBATCH --job-name=kaiju_ant
#SBATCH --output=logs/kaiju_antt%A.out
#SBATCH --error=logs/kaiju_antt%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem 16G

# Load environment (adjust if needed)
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

INPUT_DIR="results/05_kaiju"

for file in $INPUT_DIR/*.out; do
    sample=$(basename "$file" .kaiju.out)
    echo "Annotating rank $sample"

    awk '$1 == "C" {print $3}' "$file" | sort -u > "${INPUT_DIR}/${sample}_taxids.txt"

    taxonkit lineage "${INPUT_DIR}/${sample}_taxids.txt" -r -n > "${INPUT_DIR}/${sample}_taxid_lineage.txt"

    # Build lookup table with taxid -> lineage + rank
    awk -F'\t' '{print $1 "\t" $2 "\t" $4}' "${INPUT_DIR}/${sample}_taxid_lineage.txt" > "${INPUT_DIR}/${sample}_taxid_lineage_rank.tsv"

    awk 'BEGIN {
            FS = OFS = "\t";
            while (getline < "'"${INPUT_DIR}/${sample}_taxid_lineage_rank.tsv"'") {
                taxid[$1] = $2 "\t" $3;
            }
        }
        {
            if ($1 == "U" || $3 == "0") {
                print $0, "Unclassified", "Unclassified";
            } else {
                info = (taxid[$3] ? taxid[$3] : "NA\tNA");
                print $0, info;
            }
        }' "$file" > "${INPUT_DIR}/${sample}_with_lineage_and_rank.txt"

    echo "Finished annotating $sample"
done


# Loop over output files to generate kraken-style report
for file in $INPUT_DIR/*with_lineage_and_rank.txt; do
  sample=$(basename "$file" _with_lineage_and_rank.txt)
  echo "Processing $sample"
  python3 scripts/krakenize_kaiju.py "$file" "$INPUT_DIR/${sample}_kraken_report.txt"
done

# Clean temp files
cd $INPUT_DIR
rm *lineage_rank.tsv *taxid_lineage.txt *taxids.txt *with_lineage_and_rank.txt
