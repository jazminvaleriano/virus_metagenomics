#!/bin/bash
#SBATCH --job-name=kaiju
#SBATCH --output=logs/kaiju_%A.out
#SBATCH --error=logs/kaiju_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate kaiju_env

DB_DIR="databases/kaiju"
DB="$DB_DIR/nr_euk/kaiju_db_nr_euk.fmi"
READS_DIR="filtered_reads"
OUT_DIR="results/kaiju"

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq.gz; do
    sample=$(basename "$fq" _filtered.fastq.gz)

    # Run Kaiju classification
    kaiju \
        -t "$DB_DIR/nodes.dmp" \
        -f "$DB" \
        -i "$fq" \
        -z 4 \
        -o "$OUT_DIR/${sample}.kaiju.out" \
        -x

    # Generate summary table
    kaiju2table \
        -t "$DB_DIR/nodes.dmp" \
        -n "$DB_DIR/names.dmp" \
        -e \
        -r species \
        -o "$OUT_DIR/${sample}.kaiju.summary" \
        "$OUT_DIR/${sample}.kaiju.out"
done

conda deactivate