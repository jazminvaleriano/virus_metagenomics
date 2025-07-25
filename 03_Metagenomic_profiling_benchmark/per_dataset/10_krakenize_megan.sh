#!/bin/bash
#SBATCH --job-name=megan_tab-array
#SBATCH --output=logs/megan_tab_%a.out
#SBATCH --error=logs/megan_tab_%a.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=125G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch


module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

export PATH="/storage/research/vetsuisse_ivi/jvaleriano/tools:$PATH"

MEGAN_DIR=results/06_diamond-megan

PATH_TABLE=$MEGAN_DIR/path_table.tsv

# Generate path_table
rma-tabuliser -d $MEGAN_DIR -p -u -V -t ${SLURM_CPUS_PER_TASK}
mv $MEGAN_DIR/count_table.tsv $PATH_TABLE

# Clean up: keep dynamic headers and add dummy taxids
head -n 1 "$PATH_TABLE" | awk -F'\t' '{print "taxid\tpath\t" substr($0, index($0, $2))}' > "$MEGAN_DIR/header_with_path.txt"
tail -n +2 "$PATH_TABLE" | awk -F'\t' 'BEGIN {OFS="\t"} {print NR, $0}' > "$MEGAN_DIR/body_with_taxid.tsv"

cat "$MEGAN_DIR/header_with_path.txt" "$MEGAN_DIR/body_with_taxid.tsv" > "$MEGAN_DIR/final_node_with_path.tsv"

# Clean up temporary files
rm "$MEGAN_DIR/header_with_path.txt" "$MEGAN_DIR/body_with_taxid.tsv"

# Loop over all MEGAN output files
for file in $MEGAN_DIR/*.rma6; do
  sample=$(basename "$file" .rma6)
  echo "Processing $sample"
  python3 scripts/krakenize_megan.py "$sample"
done
