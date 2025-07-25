#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH --output=logs/centrifuge_%A.out
#SBATCH --error=logs/centrifuge_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=520G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

DB=/storage/research/vetsuisse_ivi/jvaleriano/Other_databases/centrifuge_nr2018/nt
READS_DIR=00_simulated_reads
OUT_DIR="results/07_centrifuge"

mkdir -p "$OUT_DIR"

for fq in $READS_DIR/*.fastq; do
    sample=$(basename "$fq" .fastq)

    start_time=$(date +%s)

    echo "[$(date)] Starting classification for $sample"

    centrifuge \
        -x $DB \
        -U $fq \
        --report-file "$OUT_DIR/${sample}.centrifuge.report.tsv" \
        -S $OUT_DIR/${sample}.centrifuge.out \
        --threads 8
    
    centrifuge-kreport -x $DB $OUT_DIR/${sample}.centrifuge.out > $OUT_DIR/${sample}.centrifuge.kreport

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    echo "[$(date)] Finished $sample in $elapsed seconds"

    # Log to CSV
    echo "kraken,$sample,$start_time,$end_time,$elapsed" >> "runtime_benchmark.csv"

done

conda deactivate