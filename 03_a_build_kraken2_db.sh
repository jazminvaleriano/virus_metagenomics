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

module load Kraken2/2.1.2-gompi-2021a

DB_DIR="databases/kraken2_custom"

# Clean any previous partial builds
rm -rf "$DB_DIR"
mkdir -p "$DB_DIR"

# Download taxonomy and specific libraries (bacteria + viral + human)
kraken2-build --download-taxonomy --db "$DB_DIR"
kraken2-build --download-library bacteria --db "$DB_DIR"
kraken2-build --download-library viral --db "$DB_DIR"
kraken2-build --download-library human --db "$DB_DIR"

# Build the database
kraken2-build \
  --build \
  --threads 32 \
  --max-db-size 450000000000 \
  --db "$DB_DIR"

# Cleanup to save space
kraken2-build --clean --db "$DB_DIR"
