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

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuger_env

# -------------------- CONFIG --------------------
CENTRIFUGER_DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/centrifuger_RefSeqViral/centrifuger_RefSeqViral
ENRICHED_READS_DIR=02_virus_enriched_reads/centrifuge
READS_DIR=00_simulated_reads 
OUT_DIR=1st_classification_results/centrifuger
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
      --cl "${OUT_DIR}/${sample}_classified_reads" \
      > ${OUT_DIR}/${sample}report.tsv


done

echo "Centrifuger run complete for all samples."

# Move enriched reads to separate folder
mkdir $ENRICHED_READS_DIR
mv $OUT_DIR/*_classified_reads.fq.gz $ENRICHED_READS_DIR

echo "Enriched reads saved in ${ENRICHED_READS_DIR}."

conda deactivate
