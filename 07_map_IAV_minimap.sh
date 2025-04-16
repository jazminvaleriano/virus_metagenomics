#!/bin/bash
#SBATCH --job-name=SIV_alignment
#SBATCH --output=logs/SIV_align%A.out
#SBATCH --error=logs/SIV_align%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

module load minimap2/2.20-GCCcore-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# Paths
#REFERENCE="references/influenza/influenza_custom.mmi" # If using custom NCBI
REFERENCE="references/influenza/swIAVsequences.mmi"

#READS_DIR="filtered_reads"  # Tried with already filtered out host reads (pig)
#READS_DIR="raw_data" # Tried with all reads
READS_DIR="reads_targeted"

OUT_DIR="results/swIAVsequences_targeted_map"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq; do #fastq.gz if compressed
    sample=$(basename "$fq" .fastq) #fastq.gz if compressed

    # Align filtered Nanopore reads to Influenza Virus reference
    minimap2 -t 4 -ax map-ont "$REFERENCE" "$fq" | \
        samtools view -@ 4 -bS - | \
        samtools sort -@ 4 -o "$OUT_DIR/${sample}_SIV.sorted.bam"

    # Index BAM file
    samtools index "$OUT_DIR/${sample}_SIV.sorted.bam"

    # Generate alignment statistics
    samtools flagstat "$OUT_DIR/${sample}_SIV.sorted.bam" > "$OUT_DIR/${sample}_SIV.stats.txt"
    samtools idxstats "$OUT_DIR/${sample}_SIV.sorted.bam" > "$OUT_DIR/${sample}_SIV.coverage.txt"
     #coverage: reference sequence name, sequence length, # mapped read-segments and # unmapped read-segments. 

done