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
PIG_GENOME=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references/pig/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
HUMAN_GENOME=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references/human/GRCh38.primary_assembly.genome.fa
READS_DIR=$WORKDIR/01_length-filtered_reads
OUT_DIR=$WORKDIR/02_host_depleted_reads

mkdir -p "$OUT_DIR"

# Index genomes (run just once)
if [ ! -f "$PIG_GENOME.mmi" ]; then
    minimap2 -d "$PIG_GENOME.mmi" "$PIG_GENOME"
fi

if [ ! -f "$HUMAN_GENOME.mmi" ]; then
    minimap2 -d "$HUMAN_GENOME.mmi" "$HUMAN_GENOME"
fi

# Process reads
for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)

    # Step 1: map to pig genome, keep unmapped
    samtools fastq -@ 4 <(minimap2 -t 4 -ax map-ont "$PIG_GENOME.mmi" "$fq" | samtools view -@ 4 -b -f 4 -) > "$OUT_DIR/${sample}_step1_unmapped.fastq"

    # Step 2: map to human genome, keep unmapped
    samtools fastq -@ 4 <(minimap2 -t 4 -ax map-ont "$HUMAN_GENOME.mmi" "$OUT_DIR/${sample}_step1_unmapped.fastq" | samtools view -@ 4 -b -f 4 -) > "$OUT_DIR/${sample}_no_host.fastq"

    # Cleanup
    rm "$OUT_DIR/${sample}_step1_unmapped.fastq"
done