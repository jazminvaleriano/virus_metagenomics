#!/bin/bash
#SBATCH --job-name=kaiju_db
#SBATCH --output=logs/kaiju_db_%A.out
#SBATCH --error=logs/kaiju_db_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate kaiju_env

# Set paths
DB_DIR="databases/kaiju_std"
TMP_DIR="$DB_DIR/tmp"
KRAKEN_DB="databases/kraken2_std"
KRAKEN_MAPFILE="$KRAKEN_DB/seqid2taxid.map"
TAX_TREE="$KRAKEN_DB/taxonomy/nodes.dmp"
TAX_NAMES="$KRAKEN_DB/taxonomy/names.dmp"

# Create directories
mkdir -p "$DB_DIR" "$TMP_DIR"
cd "$TMP_DIR"

echo "=== Downloading RefSeq proteins for Bacteria, Archaea, Viruses, Human ==="

# Download RefSeq proteins
wget -r -nd -A '*.faa.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/
wget -r -nd -A '*.faa.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/archaea/
wget -r -nd -A '*.faa.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/
wget -O human.protein.faa.gz ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/human.protein.faa.gz

echo "=== Decompressing files ==="
gunzip *.faa.gz

echo "=== Concatenating protein sequences ==="
cat *.faa > "$DB_DIR/kaiju_std.faa"

echo "=== Copying taxonomy files from Kraken DB ==="
cp "$TAX_TREE" "$DB_DIR/nodes.dmp"
cp "$TAX_NAMES" "$DB_DIR/names.dmp"

cd "$DB_DIR"

echo "=== Building Kaiju index ==="
kaiju-mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o kaiju_std kaiju_std.faa
kaiju-mkfmi kaiju_std

echo "=== Done: Kaiju standard-equivalent DB built at $DB_DIR ==="

conda deactivate
