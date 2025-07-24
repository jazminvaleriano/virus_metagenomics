#!/bin/bash
#SBATCH --job-name=map_viral_reads
#SBATCH --error=logs/map_viral_reads_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

## Remapping reads classified as viral to selected reference genomes

# Load necessary modules
module load minimap2/2.20-GCCcore-10.3.0
module load SAMtools/1.13-GCC-10.3.0
THREADS=${SLURM_CPUS_PER_TASK}

# -------- CONFIGURATION --------
GENOMES_DIR=01_classification_results/viral_hits_per_sample
READS_DIR=00_reads
OUTDIR=02_remapping_results
RENAMED_FASTQ_DIR=${OUTDIR}/renamed_fastq

mkdir -p "$OUTDIR" "$RENAMED_FASTQ_DIR"

# -------- PROCESS EACH SAMPLE --------
for RAW_FASTQ in "$READS_DIR"/*.fastq; do
    SAMPLE=$(basename "$RAW_FASTQ" .fastq)
    REF="${GENOMES_DIR}/${SAMPLE}_viral_hits.fasta"
    RENAMED_FASTQ="${RENAMED_FASTQ_DIR}/${SAMPLE}_renamed.fastq"
    BAM="${OUTDIR}/${SAMPLE}.bam"

    # Skip if reference does not exist
    if [[ ! -f "$REF" ]]; then
        echo "Reference genome not found for $SAMPLE. Skipping..."
        continue
    fi

    echo "Renaming reads for $SAMPLE..."
    awk 'NR % 4 == 1 {printf "@read_%07d\n", ++i} NR % 4 != 1 {print}' "$RAW_FASTQ" > "$RENAMED_FASTQ"

    echo "Mapping sample $SAMPLE..."
    minimap2 -t "$THREADS" -ax map-ont "$REF" "$RENAMED_FASTQ" \
        | samtools sort -@ "$THREADS" -o "$BAM"

    samtools index "$BAM"
done

echo "All samples processed."
