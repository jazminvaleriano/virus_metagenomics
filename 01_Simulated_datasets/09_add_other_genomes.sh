#!/bin/bash
#SBATCH --job-name=dl_extra_genomes
#SBATCH --output=logs/dl_extra_genomes_%A.out
#SBATCH --error=logs/dl_extra_genomes_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# I decided to add some viruses to the simulated set. 
# NOTE: The list (other_species_taxids.tsv)needs to be created manually. 

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# -------- CONFIG --------
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets
EXTRA_SPECIES_FILE=$WORKDIR/viral_training_set/other_species_taxids.tsv
GENOMES_DIR=$WORKDIR/ref_genomes
GENOME_LIST_TRAINING=$GENOMES_DIR/genome_list_training.tsv
GENOME_LIST_EXTRA=$GENOMES_DIR/genome_list_extra.tsv

OUTDIR=$WORKDIR/viral_simulation_config
GENOME_LIST_CONCAT=$OUTDIR/genome_list_concat.tsv
LOG_MISSING=$OUTDIR/missing_genomes_extra.log

mkdir -p "$OUTDIR"
> "$GENOME_LIST_EXTRA"
> "$LOG_MISSING"
# ------------------------

mapfile -t extra_species_lines < "$EXTRA_SPECIES_FILE"

for line in "${extra_species_lines[@]}"; do
    [[ -z "$line" ]] && continue
    species=$(echo "$line" | cut -f1)
    taxid=$(echo "$line" | cut -f2)

    echo "Processing: $species (TaxID: $taxid)"

    ftp=$(esearch -db assembly -query "txid${taxid}[Organism:exp]" \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_RefSeq \
        | grep -v "na" | head -n1)

    if [ -z "$ftp" ]; then
        echo "$species ($taxid) - NOT FOUND" | tee -a "$LOG_MISSING"
        continue
    fi

    file=$(basename "$ftp")
    genome_path="$GENOMES_DIR/${file}_genomic.fna.gz"
    url="${ftp}/${file}_genomic.fna.gz"

    echo "→ Downloading: $species"
    wget -q -O "$genome_path" "$url"
    
    # Decompress the genome
    gunzip -f "$genome_path"
    genome_path="${genome_path%.gz}"

    echo -e "${species}\t${genome_path}" >> "$GENOME_LIST_EXTRA"
done

# Concatenate both genome lists
cat "$GENOME_LIST_TRAINING" "$GENOME_LIST_EXTRA" > "$GENOME_LIST_CONCAT"

echo "✅ Done. Combined genome list saved to: $GENOME_LIST_CONCAT"
