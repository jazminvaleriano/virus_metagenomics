#!/bin/bash

DEPTH_DIR=03_coverage_metrics/depth_files
STATS_DIR=03_coverage_metrics/accession_mapping_stats
TAXON_MAP=01_classification_results/viral_hits_per_sample/filtered_accession2taxid.tsv
OUTDIR=04_positive_results
SCRIPT=scripts/calculate_viral_metrics.py

mkdir -p "$OUTDIR"

for DEPTH_FILE in "$DEPTH_DIR"/*_depth.txt; do
    SAMPLE=$(basename "$DEPTH_FILE" _depth.txt)

    echo "▶️ Processing $SAMPLE..."

    python "$SCRIPT" \
        "$DEPTH_FILE" \
        "$STATS_DIR/${SAMPLE}_accession_stats.tsv" \
        "$TAXON_MAP" \
        "$OUTDIR/${SAMPLE}_viral_metrics.tsv"
done

echo "✅ All samples processed."
