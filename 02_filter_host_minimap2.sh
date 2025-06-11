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
GENOME="references/pig/pig.mmi"
READS_DIR="/storage/research/vetsuisse_ivi/jvaleriano/20250324_1831_MN34349_FAY63474_1833a9d4/Rerun_basecall_20250603/fastq"
OUT_DIR="filtered_reads/second_basecall"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)
    
    # Map reads to pig genome
    minimap2 -t 4 -ax map-ont "$GENOME" "$fq" | \
    samtools view -@ 4 -b -f 4 -o "$OUT_DIR/${sample}_unmapped.bam" -

    # Convert unmapped BAM to FASTQ
    samtools fastq -@ 4 "$OUT_DIR/${sample}_unmapped.bam" | gzip > "$OUT_DIR/${sample}_filtered.fastq"

    # Clean up BAM
    rm "$OUT_DIR/${sample}_unmapped.bam"
done
