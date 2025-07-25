#!/bin/bash
#SBATCH --job-name=diamond_array
#SBATCH --output=logs/diamond_%A_%a.out
#SBATCH --error=logs/diamond_%A_%a.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=520G
#SBATCH --array=0-4   #0-<N> Where N is the number of FASTQ files minus 1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# To build reads list can run: 
#ls 00_simulated_reads/*.fastq > fastq_list.txt

# -------------------- CONFIG --------------------
READS_LIST=fastq_list.txt
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/DIAMOND_nr/nr
OUT_DIR=results/04_diamond
THREADS=${SLURM_CPUS_PER_TASK}
BLOCK_SIZE=15 # Bigger --block-size increases use of mem, but improve performance. Expect to use ~ 6x this memory (in GB), but never go >20

# ------------------------------------------------


mkdir -p "$OUT_DIR"

# ==== GET FASTQ ====
FQ=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $READS_LIST)
SAMPLE=$(basename "$FQ" .fastq)

start_time=$(date +%s)

echo "Running DIAMOND LCA for $SAMPLE"
echo "[$(date)] Starting classification for $SAMPLE"

# If using the option --fast, identity is >90%, default is 60%
diamond blastx \
    --db "$DB_DIR" \
    --query "$FQ" \
    --out "$OUT_DIR/${SAMPLE}.out" \
    --outfmt 102 \
    --include-lineage \
    --threads "$THREADS" \
    --fast \
    --block-size $BLOCK_SIZE \
    --index-chunks 4

    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    echo "[$(date)] Finished $SAMPLE in $elapsed seconds"

    # Log to CSV
    echo "diamond,$SAMPLE,$start_time,$end_time,$elapsed" >> "runtime_benchmark.csv"

conda deactivate
