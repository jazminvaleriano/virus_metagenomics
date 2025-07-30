#!/bin/bash

# Define input and output paths
MEDAKA_DIR="05_variant_calling/medaka"
REF_DIR="/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references"

# Load required modules
module load BCFtools/1.13-GCC-10.3.0
module load HTSlib/1.12-GCC-10.3.0

# Loop through each barcode-virus combination
for sample_dir in "$MEDAKA_DIR"/*/; do
    sample_virus=$(basename "$sample_dir")
    barcode=${sample_virus%%_*}
    virus=${sample_virus#*_}

    VCF="$sample_dir/medaka.sorted.vcf"
    VCF_GZ="${VCF}.gz"
    REF="$REF_DIR/$virus/$virus.fasta"
    OUT="$sample_dir/${sample_virus}_consensus.fasta"

    # Check input files
    if [[ ! -f "$VCF" ]]; then
        echo "VCF not found: $VCF"
        continue
    fi
    if [[ ! -f "$REF" ]]; then
        echo "Reference not found: $REF"
        continue
    fi

    # Compress and index VCF if needed
    if [[ ! -f "$VCF_GZ" ]]; then
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
    fi

    # Generate consensus sequence
    bcftools consensus -f "$REF" "$VCF_GZ" > "$OUT"

    echo "Consensus written to: $OUT"
done