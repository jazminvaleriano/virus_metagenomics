#!/bin/bash

# Paths 
NP_REPORTS_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/qc"
SUMMARY_FILE="$NP_REPORTS_DIR/summary_original.tsv"

# Create table header
echo -e "Barcode\tMeanReadLength\tMeanReadQuality\tMedianReadLength\tMedianReadQuality\tNumReads\tReadLengthN50\tSTDEV_ReadLength\tTotalBases\tReads>Q10\tReads>Q15\tReads>Q20\tReads>Q25\tReads>Q30" > "$SUMMARY_FILE"

# Iterate through each nanoplot report
for report in "$NP_REPORTS_DIR"/barcode*/NanoStats.txt; do
    barcode=$(basename "$(dirname "$report")")

    # Extract values with grep and awk
    mean_len=$(grep "Mean read length" "$report" | awk -F ':' '{gsub(/,/, "", $2); print $2}' | xargs)
    mean_qual=$(grep "Mean read quality" "$report" | awk -F ':' '{print $2}' | xargs)
    median_len=$(grep "Median read length" "$report" | awk -F ':' '{gsub(/,/, "", $2); print $2}' | xargs)
    median_qual=$(grep "Median read quality" "$report" | awk -F ':' '{print $2}' | xargs)
    num_reads=$(grep "Number of reads" "$report" | awk -F ':' '{print $2}' | xargs)
    n50=$(grep "Read length N50" "$report" | awk -F ':' '{gsub(/,/, "", $2); print $2}' | xargs)
    stdev_len=$(grep "STDEV read length" "$report" | awk -F ':' '{print $2}' | xargs)
    total_bases=$(grep "Total bases" "$report" | awk -F ':' '{gsub(/,/, "", $2); print $2}' | xargs)

    q10=$(grep "^>Q10" "$report" | awk '{print $2}')
    q15=$(grep "^>Q15" "$report" | awk '{print $2}')
    q20=$(grep "^>Q20" "$report" | awk '{print $2}')
    q25=$(grep "^>Q25" "$report" | awk '{print $2}')
    q30=$(grep "^>Q30" "$report" | awk '{print $2}')

    # Write line to file
    echo -e "${barcode}\t${mean_len}\t${mean_qual}\t${median_len}\t${median_qual}\t${num_reads}\t${n50}\t${stdev_len}\t${total_bases}\t${q10}\t${q15}\t${q20}\t${q25}\t${q30}" >> "$SUMMARY_FILE"
done
