#!/bin/bash
#SBATCH --job-name=DIAMOND_run
#SBATCH --output=logs/diamond_run%A.out
#SBATCH --error=logs/diamond_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load modules
module load DIAMOND/2.1.8-GCC-10.3.0

# -------------------- CONFIG --------------------
CLASS_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/results/centrifuger_prebuiltIDX/"
DB_PREFIX="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/databases/centrifuger_prebuilt/cfr_gtdb_r226+refseq_hvfc"
# ------------------------------------------------

cd "$CLASS_DIR"

echo "Starting quantification..."

for f in barcode*; do
    sample=$(basename "$f")
    echo "Processing $sample"
    
    centrifuger-quant \
      -x "$DB_PREFIX" \
      -c "$f" \
      > "${sample}_report.tsv"
done

echo "Quantification done for all samples."
conda deactivate
