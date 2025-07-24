#!/bin/bash
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch


ACCESSION2TAXID=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/taxonomy/nucl_gb.accession2taxid
OUTPUT_DIR=01_classification_results/viral_hits_per_sample
FILTERED_DB=${OUTPUT_DIR}/filtered_accession2taxid.tsv

# 1. Collect all unique taxIDs from all samples
echo "Collecting all unique taxIDs..."
cat "$OUTPUT_DIR"/*_hits_taxid.tsv | cut -f2 | sort -u > "$OUTPUT_DIR/all_taxids.tsv"

# 2. Filter the massive accession2taxid file only once
echo "Filtering accession2taxid (this may take a few minutes)..."
awk -F"\t" 'NR==FNR {wanted[$1]; next} ($3 in wanted && $1 ~ /^(NC_|NZ_)/)' "$OUTPUT_DIR/all_taxids.tsv" "$ACCESSION2TAXID" > "$FILTERED_DB"

echo "✅ Filtered accession2taxid saved to $FILTERED_DB"
