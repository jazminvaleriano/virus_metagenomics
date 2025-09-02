#!/bin/bash

FASTA_DIR="05_variant_calling/medaka"
OUTPUT="accession_description.tsv"

# Empty or create the output file
> "$OUTPUT"

# Loop through all FASTA files
find "$FASTA_DIR" -type f -name "*.fasta" | while read -r fasta; do
  grep "^>" "$fasta" | while read -r header; do
    accession=$(echo "$header" | cut -d' ' -f1 | sed 's/^>//')
    description=$(echo "$header" | cut -d' ' -f2-)
    echo -e "$accession\t$description" >> "$OUTPUT"
  done
done

echo "Dictionary saved to $OUTPUT"
