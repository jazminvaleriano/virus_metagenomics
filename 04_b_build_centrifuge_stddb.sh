#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH --output=logs/centrifuge_db%A.out
#SBATCH --error=logs/centrifuge_db%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# This script builds a centrifuge db comparable to Kraken's "standard database"
# using NCBI sequences downloaded in previous step. 

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuge_env

module load parallel/20230722-GCCcore-12.2.0
set -euo pipefail

# Paths to Kraken files
KRAKEN_DB="databases/kraken2_std"
TAX_TREE="$KRAKEN_DB/taxonomy/nodes.dmp"
TAX_NAMES="$KRAKEN_DB/taxonomy/names.dmp"

# Centrifuge files
CENTRIFUGE_DB="databases/centrifuge_std"
CONVERSION_TABLE="$CENTRIFUGE_DB/seqid2taxid.map"

# ## STEP 1: Safely concatenate the sequences

# echo "Gathering all .fna files..."
# find "$CENTRIFUGE_DB/library" -name "*.fna" > all_fna_files.txt
# echo "$(wc -l < all_fna_files.txt) .fna files found."


# echo "Splitting file list into 5000-line chunks..."
# split -l 5000 all_fna_files.txt chunks_

# echo "Creating temp directory for partial outputs..."
# mkdir -p tmp_cat_parts

# echo "Running parallel concatenation..."
# ls chunks_* | parallel -j 12 '
#     CHUNK={};
#     PART_NUM=$(echo $CHUNK | sed "s/chunks_//");
#     xargs -a $CHUNK cat > tmp_cat_parts/part_${PART_NUM}.fna
# '

# echo "Merging partial files into final input-sequences.fna..."
# cat tmp_cat_parts/part_*.fna > "$CENTRIFUGE_DB/input-sequences.fna"

# echo "Cleaning up temporary files..."
# rm -r tmp_cat_parts chunks_* all_fna_files.txt

# ## STEP 2: Append viral sequences
# echo "Starting viral sequence concatenation..."

# # Path to viral library
# VIRAL_PATH="$CENTRIFUGE_DB/library/viral"

# # Step 0: Verify initial sequence count 
# echo "Counting total sequences in final input-sequences.fna..."
# TOTAL_COUNT=$(grep -c "^>" "$CENTRIFUGE_DB/input-sequences.fna")
# echo "Total sequences before append: $TOTAL_COUNT"

# # Step 1: Decompress and concatenate into intermediate file
# echo "Concatenating all .fna.gz files from: $VIRAL_PATH"
# find "$VIRAL_PATH" -name "*.fna.gz" | sort | xargs zcat > viral-sequences.fna

# # Step 2: Count viral sequences
# echo "Counting number of sequences in viral-sequences.fna..."
# VIRAL_COUNT=$(grep -c "^>" viral-sequences.fna)
# echo "Viral sequences: $VIRAL_COUNT"

# # Step 3: Append to main input-sequences.fna
# echo "Appending viral-sequences.fna to input-sequences.fna..."
# cat viral-sequences.fna >> "$CENTRIFUGE_DB/input-sequences.fna"

# # Step 4: Verify total sequence count
# echo "Counting total sequences in final input-sequences.fna..."
# TOTAL_COUNT=$(grep -c "^>" "$CENTRIFUGE_DB/input-sequences.fna")
# echo "Total sequences after append: $TOTAL_COUNT"

# # Step 5: Clean up
# echo "Cleaning up intermediate file..."
# rm viral-sequences.fna

# echo "🎉 Done appending viral sequences."

## STEP 3: Build Centrifuge DB
centrifuge-build \
  -p 8 \
  --bmax 33554432 \
  --conversion-table "$CONVERSION_TABLE" \
  --taxonomy-tree "$TAX_TREE" \
  --name-table "$TAX_NAMES" \
  "$CENTRIFUGE_DB/input-sequences.fna" \
  "$CENTRIFUGE_DB/centrifuge_std"


conda deactivate
