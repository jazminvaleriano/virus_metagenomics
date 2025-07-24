#!/bin/bash
#SBATCH --job-name=subsample_ds
#SBATCH --output=logs/subsample_ds.out
#SBATCH --error=logs/subsample_ds.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2:00:00
#SBATCH --mem=8G

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

#=====CONFIG==================================
INPUT_FILE=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/full_reads_sim_om/reads/full_sim_om.fastq
OUTPUT_DIR=00_mock_control_reads/full_reads_sim_om/reads
SEED=100
read_counts=(5000) # Subsample sizes 
#read_counts=(200000 100000 50000 20000 10000) # Subsample sizes 
# --------------------------------------------

# Create output directory
mkdir -p $OUTPUT_DIR

# Loop over each target read count
for count in "${read_counts[@]}"; do
    echo "Generating subsample with $count reads..."
    output_file="${OUTPUT_DIR}/full_sim_om_${count}.fastq"
    seqtk sample -s$SEED $INPUT_FILE $count > $output_file
done

echo "All subsamples generated in $OUTPUT_DIR/"

conda deactivate