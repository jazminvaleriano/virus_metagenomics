#!/bin/bash
#SBATCH --job-name=chopper_filter
#SBATCH --output=logs/chopper_filter_%j.out
#SBATCH --error=logs/chopper_filter_%j.err
#SBATCH --partition=epyc2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# Set variables
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/00_Preprocessing
INPUT_DIR=$WORKDIR/00_raw_reads
OUTPUT_DIR=$WORKDIR/01_length-filtered_reads

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

for fq in $INPUT_DIR/*.fastq; do
    sample=$(basename $fq .fastq)

    # Run chopper with filters
    chopper \
      -i $fq \
      -l 200 \
      --maxlength 20000 \
      -t 4 > $OUTPUT_DIR/${sample}_length-filt.fastq
done
