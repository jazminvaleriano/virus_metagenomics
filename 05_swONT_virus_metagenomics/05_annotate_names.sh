#!/bin/bash

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# ------- CONFIG -------
CENTRIFUGER_DIR=01_classification_results/centrifuger
KRAKEN_DIR=01_classification_results/kraken
OUT_DIR_CEN=01_classification_results/centrifuger/with_names
OUT_DIR_KRA=01_classification_results/kraken/with_names

mkdir -p "$OUT_DIR_CEN" "$OUT_DIR_KRA"

# ------- Function to annotate a report -------
annotate_report () {
    local input_report="$1"
    local output_report="$2"
    local taxid_field="$3"      # 1-based column index of taxID in input
    local readid_field="$4"     # 1-based column index of readID in input
    local temp_prefix="$5"

    echo "Extracting unique taxIDs from $input_report ..."
    cut -f"$taxid_field" "$input_report" | sort -u > taxids.txt

    echo "Annotating taxIDs with scientific names..."
    taxonkit lineage taxids.txt | \
    awk -F'\t' 'BEGIN{OFS="\t"} {n=split($2,a,";"); print $1, a[n]}' > taxid_names.tsv

    echo -e "taxID\tname" | cat - taxid_names.tsv > taxid_names_with_header.tsv

    awk -F'\t' -v r="$readid_field" -v t="$taxid_field" 'BEGIN{OFS="\t"} {print $r, $t}' "$input_report" > "${temp_prefix}_core.tsv"
    echo -e "readID\ttaxID" | cat - "${temp_prefix}_core.tsv" > "${temp_prefix}_core_with_header.tsv"

    csvtk join -t -f taxID taxid_names_with_header.tsv "${temp_prefix}_core_with_header.tsv" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1, $3}' > "$output_report"

    rm taxids.txt taxid_names.tsv taxid_names_with_header.tsv "${temp_prefix}_core.tsv" "${temp_prefix}_core_with_header.tsv"
    echo "Saved: $output_report"
}


# ----------- LOOP THROUGH CENTRIFUGE REPORTS -----------
echo "🌀 Annotating all Centrifuge reports..."
for report in "$CENTRIFUGER_DIR"/*report.tsv; do
    base=$(basename "$report" .tsv)
    out_report="$OUT_DIR_CEN/${base}_withNames.tsv"
    annotate_report "$report" "$out_report" 3 1 "tmp_cen_${base}"
done

# ----------- LOOP THROUGH KRAKEN REPORTS -----------
echo "🐙 Annotating all Kraken2 reports..."
for report in "$KRAKEN_DIR"/*.kraken2.out; do
    base=$(basename "$report" .kraken2.out)
    out_report="$OUT_DIR_KRA/${base}_withNames.tsv"
    annotate_report "$report" "$out_report" 3 2 "tmp_kra_${base}"
done
