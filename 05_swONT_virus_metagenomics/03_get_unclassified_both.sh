#!/bin/bash

CENTRIFUGER_DIR=01_classification_results/centrifuger
KRAKEN_DIR=01_classification_results/kraken
OUT_DIR=01_classification_results/unclassified_both

# Create the output directory if it doesn't exist
mkdir -p "$OUT_DIR"

for cen_fq in $CENTRIFUGER_DIR/*_unclassified_reads.fq.gz ; do
    sample=$(basename "$cen_fq" _unclassified_reads.fq.gz)
    kra_fq=${KRAKEN_DIR}/${sample}_unclassified.fastq
    out_fq=${OUT_DIR}/${sample}_unclassified_both.fastq

    # Extract headers from both files
    zgrep "^@.*" "$cen_fq" | cut -d' ' -f1 > centrifuge_headers.tmp
    grep "^@.*" "$kra_fq" | cut -d' ' -f1 > kraken_headers.tmp

    # Find common headers
    comm -12 <(sort centrifuge_headers.tmp) <(sort kraken_headers.tmp) > common_headers.tmp

    # Extract full reads from Centrifuge (gzipped) file for common headers
    zcat "$cen_fq" | awk '
        BEGIN {
            while ((getline h < "common_headers.tmp") > 0) {
                headers[h] = 1
            }
        }
        /^@/ {
            header = $1
            keep = header in headers
        }
        {
            if (keep) {
                print
                line++
            } else {
                line = 0
            }
        }
        line == 4 {
            line = 0
        }
    ' > "$out_fq"
done

# Clean up temporary files
rm -f centrifuge_headers.tmp kraken_headers.tmp common_headers.tmp
