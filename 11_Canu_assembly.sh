#!/bin/bash
#SBATCH --job-name=canu_assembly
#SBATCH --output=logs/canu_assembly_%A.out
#SBATCH --error=logs/canu_assembly_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load Canu module
module load canu/2.2-GCCcore-10.3.0-Java-11

# Set paths
READS="filtered_reads"
OUTDIR="results/canu_assembly"
GENOME_SIZE="5m"  # Trying to balance for bacterial + viral contigs

mkdir -p "$OUTDIR"

for fq in "$READS"/*_filtered.fastq.gz; do
    sample=$(basename "$fq" _filtered.fastq.gz)
    mkdir -p "$OUTDIR/$sample"

    echo "Running Canu on $sample..."

    canu -p "$sample" -d "$OUTDIR/$sample" \
        genomeSize=$GENOME_SIZE \
        -nanopore-raw "$fq" \
        useGrid=false \
        maxThreads=$SLURM_CPUS_PER_TASK \
        maxMemory=128
done
