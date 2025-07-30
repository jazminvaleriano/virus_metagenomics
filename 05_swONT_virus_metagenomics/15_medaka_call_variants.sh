#!/bin/bash
#SBATCH --job-name=medaka
#SBATCH --output=logs/medaka_%A_%a.out
#SBATCH --error=logs/medaka_%A_%a.err
#SBATCH --array=0-11
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# Define paths
READ_DIR=00_reads
REF_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references
OUT_DIR=05_variant_calling/medaka
SAMPLE_LIST=05_variant_calling/sample_list.csv

mkdir -p "$OUT_DIR"

# Read the correct line from sample list
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$SAMPLE_LIST")
BARCODE=$(echo "$LINE" | cut -d',' -f1)
VIRUS=$(echo "$LINE" | cut -d',' -f2)

# Define inputs
READS=$READ_DIR/${BARCODE}*.fastq
REF="$REF_DIR/${VIRUS}/${VIRUS}.fasta"
SAMPLE_NAME="${BARCODE}_${VIRUS}"
SAMPLE_OUT="$OUT_DIR/$SAMPLE_NAME"
MODEL="r1041_e82_400bps_hac_variant_v5.0.0"

# Run medaka_variant
medaka_variant \
  -i "$READS" \
  -r "$REF" \
  -o "$SAMPLE_OUT" \
  -m "$MODEL" \
  -t 4 \
  -f
