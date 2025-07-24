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

# Goal is to create a training viral dataset for NanoSim from the real dataset. 
# First step is to classify reads with stringent criteria to find a subset of potential viral reads, and 
# a list of viruses present in the samples. 

# Load Kraken2 module
module load Kraken2/2.1.2-gompi-2021a

# -------------------- CONFIG --------------------
READS_DIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Preprocessing/02_host_depleted_reads
DB_DIR=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/kraken2_RefSeqViral
OUTDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/viral_training_set
CLASSIFICATION_DIR=$OUTDIR/classification
THREADS=${SLURM_CPUS_PER_TASK}
# ------------------------------------------------

mkdir -p $OUTDIR
mkdir -p $CLASSIFICATION_DIR

for fq in $READS_DIR/*.fastq; do
    sample=$(basename $fq _filtered.fastq)
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
    cat $CLASSIFICATION_DIR/${sample}_classified.fastq >> $OUTDIR/viral_classified_reads.fastq

    # Add found species to reference genome list
    awk -F"\t" '$4 ~ /^S/ {print $6 "\t" $5 "\t" $3}' $CLASSIFICATION_DIR/${sample}.report.txt >> $OUTDIR/species_read_counts.tsv

done

## Processing species counts file
# Eliminate trailing spaces at the beginning of line
sed -E 's/^[[:space:]]+//' $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
# Select species with read count >0
awk -F"\t" '$3 != 0' $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
# Eliminate "phage" entries
grep -iv "phage" $OUTDIR/species_read_counts.tsv > tmp && mv tmp $OUTDIR/species_read_counts.tsv
# Create a list with the unique viral species
cut -f2 $OUTDIR/species_read_counts.tsv | sort -n | uniq > $OUTDIR/unique_taxids.txt
cut -f1,2 "$OUTDIR/species_read_counts.tsv" | sort -k2,2n | uniq > "$OUTDIR/unique_species_taxids.tsv"

# Clean up unnecessary files
rm $OUTDIR/barcode*
rm -R $CLASSIFICATION_DIR

