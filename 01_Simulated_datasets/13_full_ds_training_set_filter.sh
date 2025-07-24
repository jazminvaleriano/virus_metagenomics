#!/bin/bash
#SBATCH --job-name=kraken2_filt
#SBATCH --output=logs/kraken2_filt%A.out
#SBATCH --error=logs/kraken2_filt%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Goal is to create a training dataset including diverse species for NanoSim from the real dataset. 
# First step is to classify reads with stringent criteria to find a subset of potential viral reads, and 
# a list of species present in the samples. 

# Load Kraken2 module
module load Kraken2/2.1.2-gompi-2021a

# -------------------- CONFIG --------------------
READS_DIR=/storage/research/vetsuisse_ivi/jvaleriano/00_Preprocessing/02_host_depleted_reads
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kraken2_PFP
OUTDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/full_training_set
CLASSIFICATION_DIR=$OUTDIR/classification
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

mkdir -p $OUTDIR
mkdir -p $CLASSIFICATION_DIR

for fq in $READS_DIR/*.fastq; do
    sample=$(basename $fq .fastq)
    echo "Running Kraken2 on barcode: $fq"
    kraken2 \
        --db "$DB_DIR" \
        --threads $THREADS \
        --report "$CLASSIFICATION_DIR/${sample}.report.txt" \
        --output "$CLASSIFICATION_DIR/${sample}.kraken2.out" \
        --classified-out "$CLASSIFICATION_DIR/${sample}_classified.fastq" \
        --confidence 0.4 \
        $fq
    
    # Concatenate reads in one fq file
    cat $CLASSIFICATION_DIR/${sample}_classified.fastq >> $OUTDIR/full_classified_reads.fastq

    # Add found species to reference genome list
    awk -F"\t" '$4 ~ /S/ {print $6 "\t" $5 "\t" $3}' $CLASSIFICATION_DIR/${sample}.report.txt >> $OUTDIR/species_read_counts.tsv

done

## Processing species counts file
# Eliminate trailing spaces at the beginning of line
sed -E 's/^[[:space:]]+//' $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
# Select species with read count >1
awk -F"\t" '$3 > 1' $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
# Eliminate "phage, virus and Homo" entries
grep -iv "phage" $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
grep -iv "virus" $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
grep -iv "Homo" $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
# Create a list with the unique species
cut -f2 $OUTDIR/species_read_counts.tsv | sort -n | uniq > $OUTDIR/unique_taxids.txt
cut -f1,2 "$OUTDIR/species_read_counts.tsv" | sort -k2,2n | uniq > "$OUTDIR/unique_species_taxids.tsv"

awk -F'\t' '
{
    key = $1 FS $2  # Combine species and taxon as a unique key
    count = $3
    sum[key] += count
    total += count
}
END {
    for (key in sum) {
        percentage = (sum[key] / total) * 100
        printf "%s\t%d\t%.6f\n", key, sum[key], percentage
    }
}' $OUTDIR/species_read_counts.tsv | LC_ALL=C sort -t $'\t' -k4,4nr > $OUTDIR/species_read_counts_summary.tsv

# Clean up unnecessary files
rm -R $CLASSIFICATION_DIR

