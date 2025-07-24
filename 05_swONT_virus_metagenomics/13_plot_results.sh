#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch


module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env


RESULTS_DIR=04_positive_results

for f in $RESULTS_DIR/*_viral_metrics.tsv ; do 
    python scripts/plot_results.py $f
done
