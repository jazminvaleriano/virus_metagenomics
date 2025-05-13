#!/bin/bash
#SBATCH --job-name=kraken2_build
#SBATCH --output=logs/kraken2_build_%A.out
#SBATCH --error=logs/kraken2_build_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Rebuilding kraken 2 db to use with bracken

module load Kraken2/2.1.2-gompi-2021a

DB_DIR="databases/kraken2_bracken"

mkdir -p "$DB_DIR"

# Build the database
kraken2-build \
  --standard \
  --threads 32 \
  --max-db-size 450000000000 \
  --db "$DB_DIR"

# Cleanup to save space
#kraken2-build --clean --db "$DB_DIR"
