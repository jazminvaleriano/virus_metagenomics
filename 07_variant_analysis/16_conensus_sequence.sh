#!/bin/bash
set -euo pipefail

# --- Config ---
MEDAKA_DIR="02_medaka"
REF_DIR="/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references"
MIN_DEPTH=1   # positions with depth < MIN_DEPTH will be masked to 'N'

# --- Modules ---
module load SAMtools/1.13-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load HTSlib/1.12-GCC-10.3.0

# --- Loop over sample dirs like 02_medaka/barcodeXX_virusName/ ---
for sample_dir in "$MEDAKA_DIR"/*/; do
    sample_virus=$(basename "$sample_dir")
    barcode=${sample_virus%%_*}
    virus=${sample_virus#*_}

    VCF="$sample_dir/medaka.sorted.vcf"
    VCF_GZ="${VCF}.gz"
    REF="$REF_DIR/$virus/$virus.fasta"
    OUT="$sample_dir/${sample_virus}_consensus.fasta"
    NOCOV_BED="$sample_dir/${sample_virus}.nocov.bed"

    # Try common BAM names first, else fall back to any *.bam in the directory
    if [[ -f "$sample_dir/calls_to_ref.bam" ]]; then
        BAM="$sample_dir/calls_to_ref.bam"
    elif [[ -f "$sample_dir/medaka.sorted.bam" ]]; then
        BAM="$sample_dir/medaka.sorted.bam"
    else
        # pick the first BAM if present
        BAM_CANDIDATE=$(ls $sample_dir*.bam 2>/dev/null | head -n1 || true)
        if [[ -n "${BAM_CANDIDATE:-}" && -f "$BAM_CANDIDATE" ]]; then
            BAM="$BAM_CANDIDATE"
        else
            echo "[WARN] No BAM found for $sample_virus – skipping."
            continue
        fi
    fi

    # --- Checks ---
    if [[ ! -f "$VCF" ]]; then
        echo "[WARN] VCF not found: $VCF – skipping."
        continue
    fi
    if [[ ! -f "$REF" ]]; then
        echo "[WARN] Reference not found: $REF – skipping."
        continue
    fi

    # --- Index reference (faidx) if missing ---
    if [[ ! -f "${REF}.fai" ]]; then
        samtools faidx "$REF"
    fi

    # --- Compress & index VCF if needed ---
    if [[ ! -f "$VCF_GZ" ]]; then
        bgzip -c "$VCF" > "$VCF_GZ"
        tabix -p vcf "$VCF_GZ"
    fi

    # --- Make BED of positions to mask (depth < MIN_DEPTH) ---
    # -a include all positions, -d0 no depth cap
    # Convert 1-based pos to 0-based half-open BED intervals
    samtools depth -a -d0 "$BAM" \
      | awk -v min="$MIN_DEPTH" 'BEGIN{OFS="\t"} $3<min {print $1, $2-1, $2}' > "$NOCOV_BED"

    # --- Build masked consensus ---
    # -m masks positions listed in BED to 'N'
    bcftools consensus -f "$REF" -m "$NOCOV_BED" "$VCF_GZ" > "$OUT"

    echo "[OK] Consensus written to: $OUT (masked depth < $MIN_DEPTH using $NOCOV_BED)"
done