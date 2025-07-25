#!/bin/bash
#SBATCH --job-name=kaiju
#SBATCH --output=logs/kaiju_%A.out
#SBATCH --error=logs/kaiju_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=256G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch


module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate kaiju_env

DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kaiju_refseq_nr/kaiju_db_refseq_nr.fmi
READS_DIR=00_simulated_reads
OUT_DIR=results/05_kaiju
TAXONOMY_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/taxonomy


mkdir -p "$OUT_DIR"

for fq in $READS_DIR/*.fastq; do
    sample=$(basename "$fq" .fastq)
    
    start_time=$(date +%s)
    echo "[$(date)] Starting classification for $sample"


    # Run Kaiju classification
    kaiju \
        -t $TAXONOMY_DIR/nodes.dmp \
        -f "$DB" \
        -i "$fq" \
        -z 4 \
        -o "$OUT_DIR/${sample}.kaiju.out" \
        -x #enables filtering of query sequences containing low-complexity regions using the SEG algorithm from the blast+ package

    kaiju-addTaxonNames \
        -t $TAXONOMY_DIR/nodes.dmp \
        -n $TAXONOMY_DIR/names.dmp \
        -i $OUT_DIR/${sample}.kaiju.out \
        -o $OUT_DIR/${sample}.kaiju.names.out

    # Generate summary tables
    kaiju2table \
        -t $TAXONOMY_DIR/nodes.dmp \
        -n $TAXONOMY_DIR/names.dmp \
        -r species \
        -o "$OUT_DIR/${sample}.kaiju.summary_species" \
        "$OUT_DIR/${sample}.kaiju.out"

    kaiju2table \
        -t $TAXONOMY_DIR/nodes.dmp \
        -n $TAXONOMY_DIR/names.dmp \
        -r genus \
        -o "$OUT_DIR/${sample}.kaiju.summary_genus" \
        "$OUT_DIR/${sample}.kaiju.out"

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    echo "[$(date)] Finished $sample in $elapsed seconds"

    # Log to CSV
    echo "kaiju,$sample,$start_time,$end_time,$elapsed" >> "runtime_benchmark.csv"

done

conda deactivate