#!/bin/bash
#SBATCH --job-name=centrifuger_run
#SBATCH --output=logs/centrifuger_run%A.out
#SBATCH --error=logs/centrifuger_run%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment (adjust if needed)
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuger_env


# -------------------- CONFIG --------------------
READS_DIR=00_simulated_reads
CENTRIFUGER_DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/centrifugerGTDBr226+refseq
DB_PREFIX=$CENTRIFUGER_DB/cfr_gtdb_r226+refseq_hvfc
OUT_DIR=results/03_centrifuger
THREADS=${SLURM_CPUS_PER_TASK}
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

# Create output directory
mkdir -p "$OUT_DIR"

# Run Centrifuger

for fq in $READS_DIR/*.fastq; do
    sample=$(basename "$fq" .fastq )

    start_time=$(date +%s)

    echo "Running Centrifuger on sample: $sample"
    echo "[$(date)] Starting classification for $sample"

    centrifuger \
      -x "$DB_PREFIX" \
      -u $fq \
      -t "$THREADS" \
      > $OUT_DIR/${sample}report.tsv

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    echo "[$(date)] Finished $sample in $elapsed seconds"

    # Log to CSV
    echo "centrifuger,$sample,$start_time,$end_time,$elapsed" >> "runtime_benchmark.csv"

done
echo "Centrifuger run complete for all samples."

cd $OUT_DIR
echo "Starting quantification..."

for f in *report.tsv; do
    sample=$(basename "$f" report.tsv)
    echo "Processing $sample"
    
    centrifuger-quant \
      -x "$DB_PREFIX" \
      -c "$f" \
      > "${sample}_q.tsv"

    centrifuger-kreport \
      -x $DB_PREFIX  $f >  "${sample}_kreport"
done

echo "Quantification done for all samples."

echo "Annotating kreport files with custom taxRanks..."

for kreport in *_kreport; do
    sample=$(basename "$kreport" _kreport)
    qfile="${sample}_q.tsv"
    output="${sample}_kreport_annotated"

    echo "Processing $sample"

    awk -F'\t' -v OFS='\t' '
    BEGIN {
        # Mapping of taxonomic ranks to single-letter codes
        rank_map["no rank"] = "R";
        rank_map["domain"] = "R2";
        rank_map["superkingdom"] = "R2";
        rank_map["kingdom"] = "R2";
        rank_map["phylum"] = "P";
        rank_map["class"] = "C";
        rank_map["order"] = "O";
        rank_map["family"] = "F";
        rank_map["genus"] = "G";
        rank_map["species"] = "S";
    }
    FNR==NR {
        taxrank[$2] = rank_map[$3] ? rank_map[$3] : $3;
        next;
    }
    NR==1 && $4=="U" {
        print $0; next;
    }
    {
        rank = taxrank[$5] ? taxrank[$5] : "U";
        $4 = rank;
        print $0;
    }' "$qfile" "$kreport" > "$output"

done

echo "Custom annotation done for all samples."
conda deactivate
