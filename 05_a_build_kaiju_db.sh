#!/bin/bash
#SBATCH --job-name=kaiju_db
#SBATCH --output=logs/kaiju_db_%A.out
#SBATCH --error=logs/kaiju_db_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch


module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate kaiju_env

DB_DIR="databases/kaiju"

mkdir -p "$DB_DIR"
cd "$DB_DIR"
kaiju-makedb -s nr_euk 
conda deactivate