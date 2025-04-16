#!/bin/bash

FASTQ_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/raw_data"
OUTPUT_FILE="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/qc/read_lengths_all_barcodes.tsv"

# Add header
echo -e "read_id\tread_length\tbarcode\tlog_length" > "$OUTPUT_FILE"

# Loop over all fastq.gz files
for file in "$FASTQ_DIR"/barcode*.fastq.gz; do
    barcode=$(basename "$file" | cut -d'.' -f1)  # Extract 'barcode21'

    echo "Processing $file..."

    zcat "$file" | awk -v bc="$barcode" 'NR%4==1 {
        split($1, a, "@"); id=a[2]
    } 
    NR%4==2 {
        len=length($0)
        loglen=log(len)/log(10)
        printf "%s\t%d\t%s\t%.4f\n", id, len, bc, loglen
    }' >> "$OUTPUT_FILE"
done

echo "saved all read lengths with log10 column to: $OUTPUT_FILE"
