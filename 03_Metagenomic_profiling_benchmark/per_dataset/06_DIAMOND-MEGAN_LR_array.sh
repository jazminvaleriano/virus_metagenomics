#!/bin/bash
#SBATCH --job-name=megan-array
#SBATCH --output=logs/megan_%A_%a.out
#SBATCH --error=logs/megan_%A_%a.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --array=0-6       # Adjust to the number of files minus 1
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# To build reads list can run: 
ls 00_simulated_reads/*.fastq > fastq_list.txt

# ==== CONFIG ====
READS_LIST=fastq_list.txt
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/DIAMOND_nr_NoTaxonomy/nr # DB file without extension .dmnd
MEGAN_DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/MEGAN_map2022/megan-map-Feb2022.db
OUT_DIR=results/06_diamond-megan
THREADS=${SLURM_CPUS_PER_TASK}
BLOCK_SIZE=12 # Bigger --block-size increases use of mem, but improve performance. Expect to use ~ 6x this memory (in GB), but never go >20

# ==== GET FASTQ ====
FQ=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $READS_LIST)
SAMPLE=$(basename "$FQ" .fastq)

# mkdir -p "$OUT_DIR"

# # # If using the option --fast, identity is >90%, default is 60%

# # ==== Run DIAMOND ====
    
# start_time=$(date +%s)

# diamond blastx \
#     --db $DB_DIR \
#     --query "$FQ" \
#     --out "$OUT_DIR/${SAMPLE}.daa" \
#     --outfmt 100 \
#     --threads $THREADS \
#     --block-size $BLOCK_SIZE \
#     --fast 

# # ==== Run MEGAN ====
# daa-meganizer \
#     --in "$OUT_DIR/${SAMPLE}.daa" \
#     -mdb $MEGAN_DB \
#     --longReads \
#     --threads $THREADS
#     end_time=$(date +%s)
#     elapsed=$((end_time - start_time))


# # ==== Extract taxonomic assignments ====
# daa2rma --in $OUT_DIR/${SAMPLE}.daa --out $OUT_DIR/${SAMPLE}.rma6 -mdb $MEGAN_DB

# echo "[$(date)] Finished $sample in $elapsed seconds"

# # Log to CSV
# echo "megan,$FQ,$start_time,$end_time,$elapsed" >> "runtime_benchmark.csv"

rma2info -i $OUT_DIR/${SAMPLE}.rma6 -r2c Taxonomy -o $OUT_DIR/${SAMPLE}_classification.tsv

