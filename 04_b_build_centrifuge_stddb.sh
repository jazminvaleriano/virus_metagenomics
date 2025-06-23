#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH --output=logs/centrifuge_db%A.out
#SBATCH --error=logs/centrifuge_db%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# This script builds a centrifuge db comparable to Kraken's "standard database"

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuge_env

# Paths to Kraken files
KRAKEN_DB="databases/kraken2_std"
TAX_TREE="$KRAKEN_DB/taxonomy/nodes.dmp"
TAX_NAMES="$KRAKEN_DB/taxonomy/names.dmp"

# Centrifuge files
CENTRIFUGE_DB="databases/centrifuge_std"
CONVERSION_TABLE="$CENTRIFUGE_DB/seqid2taxid_centrifuge.map"

## STEP 1: Safely concatenate the sequences

echo "Gathering all .fna files..."
find library -name "*.fna" > all_fna_files.txt

echo "Splitting file list into 5000 reads chunks..."
split -l 5000 all_fna_files.txt file_chunk_

echo "Creating temp directory for partial outputs..."
mkdir -p tmp_cat_parts

echo "Running parallel concatenation..."
module load parallel
ls file_chunk_* | parallel -j 16 'xargs -a {} cat > tmp_cat_parts/part_{#}.fna'

echo "Merging partial files into final input-sequences.fna..."
cat tmp_cat_parts/part_*.fna > input-sequences.fna

echo "Cleaning up temporary files..."
rm -r tmp_cat_parts file_chunk_* all_fna_files.txt

## STEP 2: Build Centrifuge DB
centrifuge-build \
  -p 16 \
  --conversion-table seqid2taxid.map \
  --taxonomy-tree "$TAX_TREE" \
  --name-table "$TAX_NAMES" \
  "input-sequences.fna" \
  centrifuge_std


conda deactivate
