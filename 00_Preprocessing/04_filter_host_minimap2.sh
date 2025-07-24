#!/bin/bash
#SBATCH --job-name=minimap2
#SBATCH --output=logs/minimap2%A.out
#SBATCH --error=logs/minimap2%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load minimap2/2.20-GCCcore-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# Paths
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/00_Preprocessing
GENOME=$WORKDIR/references/pig/pig.mmi
READS_DIR=$WORKDIR/01_length-filtered_reads
OUT_DIR=$WORKDIR/02_host_depleted_reads

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)
    
    # Map reads to pig genome
    minimap2 -t 4 -ax map-ont "$GENOME" "$fq" | \
    samtools view -@ 4 -b -f 4 -o "$OUT_DIR/${sample}_unmapped.bam" -

    # Convert unmapped BAM to FASTQ
    samtools fastq -@ 4 "$OUT_DIR/${sample}_unmapped.bam" > "$OUT_DIR/${sample}_filtered.fastq"

    # Clean up BAM
    rm "$OUT_DIR/${sample}_unmapped.bam"
done
