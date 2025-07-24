#!/bin/bash
#SBATCH --job-name=simulate_nanosim
#SBATCH --output=logs/simulate_nanosim_%A.out
#SBATCH --error=logs/simulate_nanosim_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load and activate environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env  

# -------- CONFIG --------
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets
GENOME_LIST=$WORKDIR/viral_simulation_config/genome_list_concat.tsv
ABUNDANCE_FILE=$WORKDIR/viral_simulation_config/abundances_simulation.tsv #Created manually
MODEL_DIR=$WORKDIR/nanosim_training_profile/nanosim_training_profile  
OUTPUT_PREFIX=$WORKDIR/00_mock_control_reads/viral_reads_sim/viral_sim
THREADS=4
# ------------------------

mkdir -p "$(dirname "$OUTPUT_PREFIX")"

# Run NanoSim simulation
simulator.py metagenome \
    -gl "$GENOME_LIST" \
    -a "$ABUNDANCE_FILE" \
    -c "$MODEL_DIR" \
    -o "$OUTPUT_PREFIX" \
    --fastq \
    -t $THREADS

echo "Simulation completed. Output: ${OUTPUT_PREFIX}.*"
