#!/bin/bash
#SBATCH --job-name=centrifuger_run
#SBATCH --output=logs/centrifuger_run%A.out
#SBATCH --error=logs/centrifuger_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment (adjust if needed)
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuger_env

# -------------------- CONFIG --------------------
READS_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/filtered_reads/second_basecall"
CENTRIFUGER_DB="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/databases/centrifuger_prebuilt/cfr_gtdb_r226+refseq_hvfc"
OUT_DIR="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/results/centrifuger"
SHEET="$OUT_DIR/sample_sheet.tsv"
THREADS=${SLURM_CPUS_PER_TASK}
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

# Create output directory
mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Create sample sheet
echo "Generating sample sheet..."
> "$SHEET"
for f in "$READS_DIR"/*.fastq; do
    base=$(basename "$f" _filtered.fastq)
    echo -e "$f . . . $base" >> "$SHEET"
done

# Run Centrifuger with full output
echo "Running centrifuger on all samples..."
centrifuger \
  -x "$CENTRIFUGER_DB" \
  --sample-sheet "$SHEET" \
  -t "$THREADS" \
  --un unclassified_reads \
  --cl classified_reads

# Output files per sample will be:
# - barcode21.tsv
# - unclassified_reads_barcode21.fq
# - classified_reads_barcode21.fq
# ...

echo "All samples classified with full output."
conda deactivate
