#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --error=logs/blast%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=250G

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

UNCLASSIFIED_DIR=01_classification_results/unclassified_both
BLAST_DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/BLAST_viralRefSeq/blast_viral_db
OUT_DIR=01_classification_results/results_unclassified_blast

mkdir -p "$OUT_DIR"

for fq in "$UNCLASSIFIED_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)

    # Convert FASTQ to FASTA on the fly and run BLAST
    seqtk seq -A "$fq" | \
    blastn -query - \
        -db "$BLAST_DB" \
        -out "${OUT_DIR}/${sample}_viral_hits.tsv" \
        -evalue 1e-20 \
        -perc_identity 90 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -max_target_seqs 3 \
        -num_threads 8
done

