#!/bin/bash

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# ------- CONFIG -------
CENTRIFUGE_REPORT=1st_classification_results/centrifuger/full_sim_om_200000report.tsv
CENTRIFUGE_OUT=1st_classification_results/centrifuger/full_sim_om_200000report_wLineage.tsv

DIAMOND_REPORT=1st_classification_results/diamond-megan/full_sim_om_200000_classification.tsv
DIAMOND_OUT=1st_classification_results/diamond-megan/full_sim_om_200000_withLineage.out

KRAKEN_REPORT=1st_classification_results/kraken/full_sim_om_200000.kraken2.out
KRAKEN_OUT=1st_classification_results/kraken/full_sim_om_200000.kraken2_withLineage.out

# -----------------------

annotate_report () {
    local input_report="$1"
    local output_report="$2"
    local taxid_field="$3"      # 1-based column index of taxID in input
    local readid_field="$4"     # 1-based column index of readID in input
    local temp_prefix="$5"

    echo "Extracting unique taxIDs from $input_report ..."
    cut -f"$taxid_field" "$input_report" | sort -u > taxids.txt

    echo " Annotating taxIDs with lineage..."
    taxonkit lineage taxids.txt | \
    taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2}' > taxid_lineage.tsv

    echo -e "taxID\tlineage" | cat - taxid_lineage.tsv > taxid_lineage_with_header.tsv

    echo " Extracting readID and taxID from report..."
    awk -F'\t' -v r="$readid_field" -v t="$taxid_field" 'BEGIN{OFS="\t"} {print $r, $t}' "$input_report" > "${temp_prefix}_core.tsv"
    echo -e "readID\ttaxID" | cat - "${temp_prefix}_core.tsv" > "${temp_prefix}_core_with_header.tsv"

    echo " Joining lineage info..."
    csvtk join -t -f taxID taxid_lineage_with_header.tsv "${temp_prefix}_core_with_header.tsv" | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $2, $1, $3}' > "$output_report"

    echo " Cleaning up..."
    rm taxids.txt taxid_lineage.tsv taxid_lineage_with_header.tsv "${temp_prefix}_core.tsv" "${temp_prefix}_core_with_header.tsv"
    echo " Saved: $output_report"
}

# ----------- RUN FOR EACH TOOL -----------
echo "CENTRIFUGE"
annotate_report "$CENTRIFUGE_REPORT" "$CENTRIFUGE_OUT" 3 1 "centrifuge"

echo "DIAMOND"
annotate_report "$DIAMOND_REPORT" "$DIAMOND_OUT" 2 1 "diamond"

echo "KRAKEN2"
annotate_report "$KRAKEN_REPORT" "$KRAKEN_OUT" 3 2 "kraken"
