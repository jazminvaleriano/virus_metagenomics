#!/bin/bash
#SBATCH --job-name=map_reads
#SBATCH --output=logs/map_reads_%A.out
#SBATCH --error=logs/map_reads_%A.err
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

GENOME_LIST=ref_genomes/genome_list_training_full.tsv
READS=full_training_set/full_classified_reads.fastq

REF_CONCATENATED_FNA=ref_genomes/all_training_genomes.fna
OUTDIR=full_training_set/mapped_reads
THREADS=${SLURM_CPUS_PER_TASK}

mkdir -p "$OUTDIR"

# 1. Combine all reference genomes
cut -f2 "$GENOME_LIST" | xargs cat > $REF_CONCATENATED_FNA
echo "All genomes concatenated into $REF_CONCATENATED_FNA"

# 2. Index the reference
minimap2 -d $OUTDIR/full_refs.mmi $REF_CONCATENATED_FNA

# 3. Map reads to viral references
minimap2 -ax map-ont -t $THREADS $OUTDIR/full_refs.mmi $READS > $OUTDIR/mapped.sam

# 4. Convert to BAM and extract mapped reads
samtools view -bS $OUTDIR/mapped.sam > $OUTDIR/mapped.bam
samtools view -b -F 4 $OUTDIR/mapped.bam > $OUTDIR/mapped_only.bam

# 5. Convert to FASTQ
samtools fastq $OUTDIR/mapped_only.bam > $OUTDIR/full_mapped_reads.fastq

# Cleanup
rm $OUTDIR/*.sam
rm $OUTDIR/*.bam


echo "Mapped reads saved to: $OUTDIR/full_mapped_reads.fastq"
