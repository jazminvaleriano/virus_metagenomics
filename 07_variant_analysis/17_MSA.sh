#!/bin/sh
#SBATCH --job-name=msa50_auto
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out

module load MAFFT/7.487-gompi-2021a-with-extensions
mafft --auto --reorder seq_alignment/seg4_bc36.fasta > seq_alignment/seg4_bc36.alignment.fasta