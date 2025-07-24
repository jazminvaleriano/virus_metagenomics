#!/bin/bash

# Run this for each set of reads changing file name

FASTQ_DIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/full_reads_sim_pretrained
OUT_DIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/read_lengths
OUTPUT_FILE=$OUT_DIR/Full_sim_pretrained.tsv

mkdir -p $OUT_DIR

# Add header
echo -e "barcode\tread_length\tlog_length" > "$OUTPUT_FILE"

# Loop over all fastq.gz files
for file in $FASTQ_DIR/*dataset_200000.fastq; do
    barcode=$(basename "$file" | cut -d'.' -f1) 

    echo "Processing $file..."

    cat "$file" | awk -v bc="$barcode" 'NR%4==2 {
        len=length($0)
        loglen=log(len)/log(10)
        printf "%s\t%d\t%.4f\n", bc, len, loglen
    }' >> "$OUTPUT_FILE"
done

echo "Saved all read lengths with log10 column to: $OUTPUT_FILE"
