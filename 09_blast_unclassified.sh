#!/bin/bash
#SBATCH --job-name=blast_unclassified
#SBATCH --output=logs/blast_unclassified_%A.out
#SBATCH --error=logs/blast_unclassified_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Activate conda environment with BLAST
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate blast_env

READS_DIR="results/kraken_unclassified_blast"

## Directories to blast vs Virus db:
#BLAST_DB="databases/blast/nt_viruses/nt_viruses" 
#OUT_DIR="results/blast/unclassified_vs_virus"

# Directories to blast vs prokaryotes:
BLAST_DB="databases/blast/nt_viruses/ref_prok" 
OUT_DIR="results/blast/unclassified_vs_prok"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*_unclassified.fastq; do
    sample=$(basename "$fq" _unclassified.fastq)

    # Convert FASTQ to FASTA (inline using awk, only the first time)
    #awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' "$fq" > "$READS_DIR/${sample}_unclassified.fasta"

    # Run BLAST on converted FASTA
    blastn \
        -query "$READS_DIR/${sample}_unclassified.fasta" \
        -db "$BLAST_DB" \
        -out "$OUT_DIR/${sample}_blast.txt" \
        -evalue 1e-5 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -num_threads 8 \
        -max_target_seqs 5

done