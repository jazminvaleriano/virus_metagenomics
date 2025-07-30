#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# Define folders
READ_DIR=00_reads     
REF_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references
OUT_DIR=05_variant_calling/minimap

mkdir -p "$OUT_DIR"

# Init sample list CSV (no header)
SAMPLE_LIST="$OUT_DIR/sample_list.csv"
> "$SAMPLE_LIST"

# Function to perform mapping
map_reads() {
  local barcode=$1
  local virus=$2

  READS="$READ_DIR/${barcode}*.fastq"
  REF="$REF_DIR/${virus}/${virus}.fasta"
  REF_INDEX="${REF%.fasta}.mmi"
  OUT_BAM="$OUT_DIR/${barcode}_${virus}.sorted.bam"

  echo "Indexing $virus ref genome"
  if [ ! -f "$REF_INDEX" ]; then
    minimap2 -d "$REF_INDEX" "$REF"
  fi
  echo "Mapping $barcode to $virus..."
  minimap2 -ax map-ont "$REF_INDEX" $READS \
    | samtools sort -@ 4 -o "$OUT_BAM"

  samtools index "$OUT_BAM"
}

# Declare barcode-virus pairs as grouped arrays
declare -A barcode_groups

barcode_groups["barcode25"]="Bocavirus Porcine_Hokovirus Porcine_partetravirus"
barcode_groups["barcode26"]="Bocavirus Porcine_Hokovirus Porcine_partetravirus"
barcode_groups["barcode33"]="Bovine_coronavirus Human_coronavirus_OC43"
barcode_groups["barcode34"]="Bovine_coronavirus Human_coronavirus_OC43"
barcode_groups["barcode35"]="influenza"
barcode_groups["barcode36"]="influenza"

# Iterate over groups
for barcode in "${!barcode_groups[@]}"; do
  viruses=(${barcode_groups[$barcode]})
  for virus in "${viruses[@]}"; do
    map_reads "$barcode" "$virus"
    echo "${barcode},${virus}" >> "$SAMPLE_LIST"
  done
done
