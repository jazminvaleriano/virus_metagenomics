#!/bin/bash

BLAST_DIR=01_classification_results/results_unclassified_blast
CENTRIFUGER_DIR=01_classification_results/centrifuger/with_names
KRAKEN_DIR=01_classification_results/kraken/with_names
REFSEQ_FNA=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/refseq_viral/viral.1.1.genomic.fna
ACCESSION2TAXID=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/taxonomy/nucl_gb.accession2taxid

OUTPUT_DIR=01_classification_results/viral_hits_per_sample
mkdir -p "$OUTPUT_DIR"

FILTER_PATTERN='phage|viruses|virinae|viricetes|virales|viridae|root|unclassified'

echo "Starting per-sample processing..."

for BLAST_FILE in "$BLAST_DIR"/*_unclassified_both_viral_hits.tsv; do
    SAMPLE=$(basename "$BLAST_FILE" _unclassified_both_viral_hits.tsv)
    echo "Processing sample: $SAMPLE"

    KRAKEN_FILE="$KRAKEN_DIR/${SAMPLE}_withNames.tsv"
    CENTRIFUGE_FILE="$CENTRIFUGER_DIR/${SAMPLE}report_withNames.tsv"

    # Skip if Kraken or Centrifuge file is missing
    if [[ ! -f "$KRAKEN_FILE" || ! -f "$CENTRIFUGE_FILE" ]]; then
        echo "Warning: Missing Kraken or Centrifuge file for $SAMPLE, skipping..."
        continue
    fi

    TMP_TAXIDS=$(mktemp)

    # Step 1: BLAST – column 13 (taxon name) → lookup in Kraken for taxID
    tail -n +2 "$BLAST_FILE" | cut -f13 | grep -Eiv "$FILTER_PATTERN" | sort -u | while read NAME; do
        TAXID=$(awk -F"\t" -v name="$NAME" 'tolower($2) == tolower(name) && $1 != "0"' "$KRAKEN_FILE" | cut -f1 | head -n 1)
        [[ -n "$TAXID" ]] && echo -e "$TAXID\t$NAME" >> "$TMP_TAXIDS"
    done

    # Step 2: Kraken – taxID (col 1), name (col 2)
    tail -n +2 "$KRAKEN_FILE" | grep -Eiv "$FILTER_PATTERN" | awk -F"\t" '$1 != "0"' | cut -f1,2 >> "$TMP_TAXIDS"

    # Step 3: Centrifuge – taxID (col 1), name (col 2)
    tail -n +2 "$CENTRIFUGE_FILE" | grep -Eiv "$FILTER_PATTERN" | awk -F"\t" '$1 != "0"' | cut -f1,2 >> "$TMP_TAXIDS"

    # Deduplicate
    HITS_LIST="${OUTPUT_DIR}/${SAMPLE}_hits_taxid.tsv"
    sort -u "$TMP_TAXIDS" > "$HITS_LIST"
    rm "$TMP_TAXIDS"
done

echo "All samples processed."
