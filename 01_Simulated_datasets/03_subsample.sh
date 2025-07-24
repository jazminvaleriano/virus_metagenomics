#!/bin/bash
#SBATCH --job-name=00_random_subsample_reads
#SBATCH --output=logs/00_random_subsample_reads.out
#SBATCH --error=logs/00_random_subsample_reads.err
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
INPUT_FILE=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/00_mock_control_reads/Zymo/SHORT/original/02_Zymo-200K_SHORT.fastq
SEED=100
read_counts=(5000 10000 20000 100000) # Subsample sizes 
# --------------------------------------------

OUTPUT_DIR=$(dirname "$INPUT_FILE")

# Loop over each target read count
for count in "${read_counts[@]}"; do
    file_name=$(basename "$INPUT_FILE" .fastq)
    echo "Generating subsample with $count reads..."
    output_file="${OUTPUT_DIR}/${file_name}_${count}.fastq"
    seqtk sample -s$SEED $INPUT_FILE $count > $output_file
done

echo "All subsamples generated in $OUTPUT_DIR/"

conda deactivate