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
#conda activate metagenomics_env  # If using my own model
conda activate NanoSim_pretrained_env # If using NanoSim's pretrained models, as their scikit version is older (0.22.1)

# -------- CONFIG --------
model_sufx=pt  # pt for pretrained, om for own model
OUTPUT_DIR=full_reads_sim_${model_sufx}
OUTPUT_PREFIX=${OUTPUT_DIR}/full_sim_${model_sufx}

MODEL_DIR=metagenome_ERR3152364_Even_v3.2.2/training
GENOME_LIST=ref_genomes/genome_list_training_full.tsv
ABUNDANCE_DESIGN=full_training_set/mock_community_proportions.tsv
TOTAL_SIZE=300000

FINAL_OUTFILE_NAME=${OUTPUT_DIR}/reads/full_sim_${model_sufx}.fastq
ABUNDANCE_FILE=$(dirname "$ABUNDANCE_DESIGN")/abundances_simulation.tsv
# ------------------------
mkdir -p "$(dirname "$OUTPUT_PREFIX")/sim_reads"

## Create abundance file
# Write header line
echo -e "Size\t$TOTAL_SIZE" > "$ABUNDANCE_FILE"

# Process each line after the header
awk -F'\t' 'NR > 1 { printf "%s\t%s\n", $1, $3 }' "$ABUNDANCE_DESIGN" >> "$ABUNDANCE_FILE"

# Run NanoSim simulation
simulator.py metagenome \
    -gl "$GENOME_LIST" \
    -a "$ABUNDANCE_FILE" \
    -c "$MODEL_DIR" \
    -o "$OUTPUT_PREFIX" \
    --fastq \
    -t $SLURM_CPUS_PER_TASK \
    --chimeric # NanoSim failed when this option was active, unless I added a dummy shrinkage rate (0.8 from mock ds)

echo "Simulation completed. Output: ${OUTPUT_PREFIX}.*"

## Rename output file and simplify read names using awk (NanoSim gives very long names to reads, which causes some tools to crash, e.g. minimap)
RAW_FASTQ=$(dirname $OUTPUT_PREFIX)/*aligned_reads.fastq

awk 'NR % 4 == 1 {printf "@read_%07d\n", ++i} NR % 4 != 1 {print}' "$RAW_FASTQ" > "$FINAL_OUTFILE_NAME"

