#!/bin/bash
#SBATCH --job-name=kraken2_run
#SBATCH --error=logs/kraken2_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# -------------------- CONFIG --------------------
READS_DIR=00_reads
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kraken2_PFP
OUT_DIR=01_classification_results/kraken
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

mkdir -p "$OUT_DIR"

for fq in "$READS_DIR"/*.fastq; do
    sample=$(basename "$fq" .fastq)

    echo "Running Kraken2 on sample: $sample"

    # kraken2 \
    #     --db "$DB_DIR" \
    #     --threads $THREADS \
    #     --report "$OUT_DIR/${sample}.report.txt" \
    #     --report-minimizer-data \
    #     --output "$OUT_DIR/${sample}.kraken2.out" \
    #     --unclassified-out "$OUT_DIR/${sample}_unclassified.fastq" \
    #     "$fq"

    echo "Annotating lineages for $sample..."

     # Extract taxids and annotate with TaxonKit
    cut -f7 "$OUT_DIR/${sample}.report.txt" | taxonkit lineage --threads $THREADS \
    | awk 'BEGIN{OFS="\t"} {print $1, $3}' \
    > "$OUT_DIR/${sample}_lineages.tsv"

    # Join lineage info back to Kraken2 report (taxid is column 7)
    awk 'BEGIN{OFS="\t"} NR==FNR {a[$1]=$2; next} {print $0, (a[$7] ? a[$7] : "NA")}' \
        "$OUT_DIR/${sample}_lineages.tsv" "$OUT_DIR/${sample}.report.txt" \
        > "$OUT_DIR/${sample}.report_with_lineage.txt"

    # Cleanup intermediate file
    rm "$OUT_DIR/${sample}_lineages.tsv"

done

echo "Kraken2 run complete."