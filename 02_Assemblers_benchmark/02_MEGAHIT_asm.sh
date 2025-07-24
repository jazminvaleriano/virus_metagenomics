#!/bin/bash
#SBATCH --job-name=megahit_asm
#SBATCH --output=logs/megahit_asm_%A_%a.out
#SBATCH --error=logs/megahit_asm_%A_%a.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load MEGAHIT/1.2.9-GCCcore-10.3.0

# Set datasets array
DATASETS_DIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/full_reads_sim_om/reads
DATASETS=( "$DATASETS_DIR"/*0.fastq )

# Select dataset for this array task
fq=${DATASETS[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "$fq" .fastq)

OUTDIR=results/02_megahit_assemblies
mkdir -p $OUTDIR

echo "Assembling $sample from $fq..."

megahit --presets meta-large --kmin-1pass -r "$fq" -o "$OUTDIR/$sample"

# Rename output file if it exists
mv $OUT_DIR/$sample/final.contigs.fa "results/02_megahit_assemblies/${sample}.fasta"

# Options suggested by megahit: 
# --kmin-1pass: if sequencing depth is low and too much memory used when build the graph of k_min
# --presets meta-large: if the metagenome is complex (i.e., bio-diversity is high, for example soil