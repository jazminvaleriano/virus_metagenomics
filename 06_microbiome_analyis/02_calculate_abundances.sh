#!/bin/bash

# Config
KRAKEN_DIR=01_classification_results/kraken
LEVEL="G"
THRESH=20 # Min number of distinct minimizers assigned to a taxon to include in microbiome analysis
OUT_DIR=02_diversity_analysis
BETA_OUT=$OUT_DIR/beta_diversity_matrix.txt

mkdir -p $OUT_DIR

# Calculate abundance based on k-mers

for f in $KRAKEN_DIR/*report_with_lineage.txt ; do
    python scripts/kmer_abundance.py $f $LEVEL $THRESH
done

# Calculate Beta-diversity

FILES_LIST=$(ls $KRAKEN_DIR/*abundance.tsv)

python scripts/beta_diversity.py -i $FILES_LIST -o $BETA_OUT




