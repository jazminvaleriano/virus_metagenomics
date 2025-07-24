#!/bin/bash
#SBATCH --job-name=centrifuger_run
#SBATCH --error=logs/centrifuger_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuger_env

# -------------------- CONFIG --------------------
CENTRIFUGER_DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/centrifuger_RefSeqViral/centrifuger_RefSeqViral
READS_DIR=00_reads 
OUT_DIR=01_classification_results/centrifuger
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------
echo "Running Centrifuger with:"
echo "READS_DIR: $READS_DIR"

# Create output directory
mkdir -p "$OUT_DIR"

# Run Centrifuger

for fq in $READS_DIR/*.fastq ; do
    sample=$(basename "$fq" .fastq)

    echo "Running Centrifuger on sample: $sample"

    centrifuger \
      -x "$CENTRIFUGER_DB" \
      -u $fq \
      -t "$THREADS" \
      --un "${OUT_DIR}/${sample}_unclassified_reads" \
      > ${OUT_DIR}/${sample}report.tsv


done

echo "Centrifuger run complete for all samples."

conda deactivate
