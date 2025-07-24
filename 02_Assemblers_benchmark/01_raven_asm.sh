#!/bin/bash
#SBATCH --job-name=raven_assembly
#SBATCH --output=logs/raven_assembly_%A_%a.out
#SBATCH --error=logs/raven_assembly_%A_%a.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# Set datasets array
DATASETS_DIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/full_reads_sim_om/reads
DATASETS=( "$DATASETS_DIR"/*0.fastq )

# Select dataset for this array task
fq=${DATASETS[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "$fq" .fastq)

OUTDIR=results/01_raven_assemblies
mkdir -p "$OUTDIR"

echo "Assembling ${sample}..."

# Run Raven with 2 polishing rounds
raven -t $SLURM_CPUS_PER_TASK -p 2 "$fq" > "$OUTDIR/${sample}.fasta"

conda deactivate
