#!/bin/bash
#SBATCH --job-name=quast
#SBATCH --output=logs/quast_%A.out
#SBATCH --error=logs/quast_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate conda environment with MetaQUAST & BLAST
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

#=========config==============
TOOL_LIST="flye raven megahit"
REF_GENOMES=/storage/research/vetsuisse_ivi/jvaleriano/02_Assemblers_evaluation/02_OM_VS_PT/simulated_DS_references

for tool in $TOOL_LIST; do

    OUT_DIR=results/03_assemblies_evaluation/$tool
    if [[ "$tool" == "raven" ]]; then
        ASSEMBLIES_DIR=results/01_raven_assemblies
    elif [[ "$tool" == "flye" ]]; then
        ASSEMBLIES_DIR=results/00_flye_assemblies
    elif [[ "$tool" == "megahit" ]]; then
        ASSEMBLIES_DIR=results/02_megahit_assemblies
    fi

    mkdir -p "$OUT_DIR"

    # Run MetaQUAST on all assemblies at once
    metaquast.py $ASSEMBLIES_DIR/*.fasta \
        -r "$REF_GENOMES" \
        -o "$OUT_DIR" \
        --threads 8

done
