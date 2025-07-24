#!/bin/bash
#SBATCH --job-name=map_viral_reads
#SBATCH --output=logs/map_viral_reads_%A.out
#SBATCH --error=logs/map_viral_reads_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

## Remapping reads classified as viral to selected reference genomes

# Load necessary modules
module load minimap2/2.20-GCCcore-10.3.0
module load SAMtools/1.13-GCC-10.3.0

# -------- CONFIGURATION --------
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets
GENOMES_DIR=$WORKDIR/ref_genomes
READS=$WORKDIR/viral_training_set/viral_classified_reads.fastq
OUTDIR=$WORKDIR/viral_training_set/mapped_reads
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "$OUTDIR"

# 1. Combine all reference genomes
cat $GENOMES_DIR/*.fna > $OUTDIR/viral_refs_combined.fasta

# 2. Index the reference
minimap2 -d $OUTDIR/viral_refs.mmi $OUTDIR/viral_refs_combined.fasta

# 3. Map reads to viral references
minimap2 -ax map-ont -t $THREADS $OUTDIR/viral_refs.mmi $READS > $OUTDIR/mapped.sam

# 4. Convert to BAM and extract mapped reads
samtools view -bS $OUTDIR/mapped.sam > $OUTDIR/mapped.bam
samtools view -b -F 4 $OUTDIR/mapped.bam > $OUTDIR/mapped_only.bam

# 5. Convert to FASTQ
samtools fastq $OUTDIR/mapped_only.bam > $OUTDIR/viral_mapped_reads.fastq

# Cleanup
rm $OUTDIR/mapped.sam

echo "Mapped reads saved to: $OUTDIR/viral_mapped_reads.fastq"
