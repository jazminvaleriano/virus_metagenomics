#!/bin/bash
#SBATCH --job-name=bracken_all
#SBATCH --output=logs/bracken_all%A.out
#SBATCH --error=logs/bracken_all%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# If omitted, the -l flag default is Species 

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate bracken_env

# Set paths
REPORT_DIR="results/kraken_vs_standard/first_basecall"
OUT_DIR="results/bracken/bracken_vs_standard_S_original" # Change if rerunning
DB="databases/kraken2_std"
READ_LEN=800
LEVEL="S"
THRESHOLD=1 # Min reads for classification

# Create output directory if not exists
mkdir -p "$OUT_DIR"

# Loop through all report files
for REPORT in "$REPORT_DIR"/*.report.txt; do
    SAMPLE=$(basename "$REPORT" .report.txt)
    echo "Running Bracken on $SAMPLE"
    bracken -d "$DB" \
            -i "$REPORT" \
            -o "$OUT_DIR/${SAMPLE}_bracken_${LEVEL}.txt" \
            -r "$READ_LEN" \
            -l "$LEVEL" \
            -t "$THRESHOLD"
done
