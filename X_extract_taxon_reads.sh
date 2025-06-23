#!/bin/bash

# Usage: scripts/X_extract_taxon_reads.sh "/path/to/reports" "taxon name" "output.csv"

REPORT_DIR="$1"
TAXON_NAME="$2"
OUTPUT_FILE="$3"

if [[ -z "$REPORT_DIR" || -z "$TAXON_NAME" || -z "$OUTPUT_FILE" ]]; then
  echo "Usage: $0 <report_directory> <taxon_name> <output_file>"
  exit 1
fi

# Ensure output file starts clean
echo "sample,taxon_name,reads" > "$OUTPUT_FILE"

# Loop through report files
for FILE in "$REPORT_DIR"/*.report.txt; do
  if [[ -f "$FILE" ]]; then
    SAMPLE=$(basename "$FILE" | cut -d'.' -f1)

    # Extract matching line(s) and get the reads column (3rd) and taxon name (last)
    MATCH=$(awk -v taxon="$TAXON_NAME" '$0 ~ taxon {print}' "$FILE")

    while IFS= read -r LINE; do
      if [[ -n "$LINE" ]]; then
        READS=$(echo "$LINE" | awk '{print $2}')
        TAXON=$(echo "$LINE" | awk '{for (i=6; i<=NF; i++) printf $i" "; print ""}' | sed 's/ *$//')
        if [[ "$TAXON" == "$TAXON_NAME" ]]; then
          echo "$SAMPLE,$TAXON,$READS" >> "$OUTPUT_FILE"
        fi
      fi
    done <<< "$MATCH"
  fi
done

echo "Done. Results written to $OUTPUT_FILE"
