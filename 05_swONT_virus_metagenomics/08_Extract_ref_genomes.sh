#!/bin/bash

# --- Paths ---
BLAST_DIR=01_classification_results/results_unclassified_blast
CENTRIFUGER_DIR=01_classification_results/centrifuger/with_names
KRAKEN_DIR=01_classification_results/kraken/with_names
REFSEQ_FNA=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/refseq_viral/viral.1.1.genomic.fna
FILTERED_DB=01_classification_results/viral_hits_per_sample/filtered_accession2taxid.tsv

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
        echo "⚠️  Warning: Missing Kraken or Centrifuge file for $SAMPLE. Skipping..."
        continue
    fi

    TMP_TAXIDS=$(mktemp)

    # --- Step 1: Extract accessions from BLAST column 2, clean and look up taxID ---
    tail -n +2 "$BLAST_FILE" \
        | cut -f2 \
        | sed -E 's/^ref\|//; s/\|$//' \
        | grep -Ev '^$' \
        | sort -u \
        | while read ACC; do
            TAXID=$(awk -v acc="$ACC" -F"\t" '$1 == acc {print $3; exit}' "$FILTERED_DB")
            [[ -n "$TAXID" && "$TAXID" -ne 0 ]] && echo -e "$TAXID\t$ACC" >> "$TMP_TAXIDS"
        done

    # --- Step 2: Add Kraken + Centrifuge taxID:name pairs (skip taxID 0) ---
    tail -n +2 "$KRAKEN_FILE" | grep -Eiv "$FILTER_PATTERN" | awk -F'\t' '$1 != 0 {print $1 "\t" $2}' >> "$TMP_TAXIDS"
    tail -n +2 "$CENTRIFUGE_FILE" | grep -Eiv "$FILTER_PATTERN" | awk -F'\t' '$1 != 0 {print $1 "\t" $2}' >> "$TMP_TAXIDS"

    # --- Step 3: Deduplicate taxIDs ---
    HITS_LIST="${OUTPUT_DIR}/${SAMPLE}_hits_taxid.tsv"
    sort -u "$TMP_TAXIDS" > "$HITS_LIST"
    rm "$TMP_TAXIDS"

    # --- Step 4: Map taxIDs to accessions ---
    GENOMES_LIST="${OUTPUT_DIR}/${SAMPLE}_hit_genomes_list.tsv"
    awk -F"\t" 'NR==FNR {wanted[$2]; next} ($3 in wanted)' "$HITS_LIST" "$FILTERED_DB" > "$GENOMES_LIST"

    # --- Step 5: Extract FASTA ---
    OUT_FASTA="${OUTPUT_DIR}/${SAMPLE}_viral_hits.fasta"
    python scripts/extract_matching_fasta.py "$GENOMES_LIST" "$REFSEQ_FNA" "$OUT_FASTA"

    echo "✅ Finished $SAMPLE:"
    echo "   - Hits list: $HITS_LIST"
    echo "   - Genomes list: $GENOMES_LIST"
    echo "   - Extracted FASTA: $OUT_FASTA"
done

echo "All samples processed."
