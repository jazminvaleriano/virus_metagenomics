#!/bin/bash
#SBATCH --job-name=flye_assembly
#SBATCH --output=logs/flye_assembly_%A_%a.out
#SBATCH --error=logs/flye_assembly_%A_%a.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-4
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load Flye
module load Flye/2.9-GCC-10.3.0

# Set datasets array
DATASETS_DIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/full_reads_sim_om/reads
DATASETS=( "$DATASETS_DIR"/*0.fastq )

# Select dataset for this array task
fq=${DATASETS[$SLURM_ARRAY_TASK_ID]}
sample=$(basename "$fq" .fastq)

OUTDIR=results/00_flye_assemblies/$sample

echo "Assembling ${sample} from ${fq}..."
mkdir -p "$OUTDIR"

# Run Flye in metagenome mode
flye --nano-hq "$fq" \
    --out-dir "$OUTDIR" \
    --meta \
    --threads "$SLURM_CPUS_PER_TASK"

# Rename output file if it exists
if [[ -f "$OUTDIR/assembly.fasta" ]]; then
    mv "$OUTDIR/assembly.fasta" "results/00_flye_assemblies/${sample}.fasta"
else
    echo "Warning: $OUTDIR/assembly.fasta not found"
fi