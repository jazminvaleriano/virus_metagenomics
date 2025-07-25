#!/bin/bash
#SBATCH --job-name=kraken_report
#SBATCH --output=logs/kraken_report_%j.out
#SBATCH --error=logs/kraken_report_%j.err
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1


# Activate environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env


RESULTS=results/04_diamond

# Loop over all DIAMOND output files
for file in $RESULTS/*.out; do
  sample=$(basename "$file" .out)
  echo "Processing $sample"
  python3 scripts/krakenize_diamond.py "$file" "$RESULTS/${sample}_kraken_report.txt"
done
