#!/bin/bash
#SBATCH --job-name=dl_genomes
#SBATCH --output=logs/dl_genomes_%A.out
#SBATCH --error=logs/dl_genomes_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# -------- CONFIG --------
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets
SPECIES_TAXID_FILE=$WORKDIR/full_training_set/mock_community_proportions.tsv
OUTDIR=$WORKDIR/ref_genomes_full
LOG_MISSING=$OUTDIR/missing_genomes_full.log
GENOME_LIST=$OUTDIR/genome_list_training_full.tsv
VIRUS_GENOME_LIST=$OUTDIR/genome_list_training.tsv

mkdir -p "$OUTDIR"
> "$GENOME_LIST"
> "$LOG_MISSING"
# ------------------------

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# Map file into array, skip header
mapfile -t species_taxid_lines < <(tail -n +2 "$SPECIES_TAXID_FILE")

for line in "${species_taxid_lines[@]}"; do
    [[ -z "$line" ]] && continue  # Skip empty lines

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
    genome_path="$OUTDIR/${file}_genomic.fna.gz"
    url="${ftp}/${file}_genomic.fna.gz"

    echo "→ Downloading: $species"
    wget -q -O "$genome_path" "$url"

    gunzip -f "$genome_path"
    genome_path="${genome_path%.gz}"

    echo -e "${species}\t${genome_path}" >> "$GENOME_LIST"
done

echo "Finished downloading genomes and building genome list."

cat $VIRUS_GENOME_LIST >> $GENOME_LIST

echo "Genome list saved at: $GENOME_LIST"
