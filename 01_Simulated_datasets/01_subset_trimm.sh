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
conda activate kaiju_env

#=====CONFIG==================================
SUBSET_NAME=SHORT
SHORT_LIM=200
LONG_LIM=20000
#======PATHS==================================
INPUT_FILE="00_mock_control_reads/Zymo/Zymo-GridION-EVEN-3Peaks-R103-merged.fastq"
TRIMMED_FILE="00_mock_control_reads/Zymo/$SUBSET_NAME/01_Zymo-trimmed_${SUBSET_NAME}.fastq"
SUBSAMPLE_FILE="00_mock_control_reads/Zymo/$SUBSET_NAME/original/02_Zymo-200K_${SUBSET_NAME}.fastq"
#=============================================

mkdir -p "00_mock_control_reads/Zymo/$SUBSET_NAME/original"

echo "Filtering reads shorter than $SHORT_LIM bp and longer than $LONG_LIM bp... as per Portik paper "
seqkit seq -m $SHORT_LIM -M $LONG_LIM $INPUT_FILE > $TRIMMED_FILE

echo "Subsampling 200K reads"
seqtk sample -s100 $TRIMMED_FILE 200000 > $SUBSAMPLE_FILE

