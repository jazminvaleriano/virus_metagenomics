#!/bin/bash
#SBATCH --job-name=nanosim_train
#SBATCH --output=logs/nanosim_train_%A.out
#SBATCH --error=logs/nanosim_train_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# -------- LOAD ENVIRONMENT --------
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# -------- CONFIG --------
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets
GENOMES_DIR=$WORKDIR/ref_genomes
GENOMES_LIST=$GENOMES_DIR/genome_list_training.tsv

READS=$WORKDIR/viral_training_set/mapped_reads/viral_mapped_reads.fastq
OUTDIR=$WORKDIR/nanosim_training_profile/nanosim_training_profile
THREADS=${SLURM_CPUS_PER_TASK}

# -------- RUN NANOSIM TRAINING --------
read_analysis.py metagenome \
  -i "$READS" \
  -gl "$GENOMES_LIST" \
  -o "$OUTDIR" \
  --fastq \
  -t "$THREADS"

echo "✅ NanoSim training completed. Output: $OUTDIR"
