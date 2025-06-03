#!/bin/bash
#SBATCH --job-name=dorado_gpu
#SBATCH --partition=gpu
#SBATCH --qos=job_gpu
#SBATCH --gpus=rtx3090:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=logs/dorado_gpu_%A.out
#SBATCH --error=logs/dorado_gpu_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load CUDA/12.6.0 

# Adding Dorado to PATH
export PATH=/storage/research/vetsuisse_ivi/jvaleriano/tools/dorado/dorado-1.0.0-linux-x64/bin:$PATH

# Paths
POD5_DIR="/storage/research/vetsuisse_ivi/jvaleriano/20250324_1831_MN34349_FAY63474_1833a9d4/pod5"
OUT_DIR="/storage/research/vetsuisse_ivi/jvaleriano/20250324_1831_MN34349_FAY63474_1833a9d4/Rerun_basecall_20250603"
MODEL_DIR="/storage/research/vetsuisse_ivi/jvaleriano/tools/dorado/models/dna_r10.4.1_e8.2_400bps_hac@v5.2.0"

mkdir -p "$OUT_DIR"

# Execute Dorado basecaller
dorado basecaller \
  --kit-name SQK-NBD114-96 \
  --min-qscore 9 \
  --emit-fastq \
  --output-dir "$OUT_DIR" \
  "$MODEL_DIR" "$POD5_DIR"
