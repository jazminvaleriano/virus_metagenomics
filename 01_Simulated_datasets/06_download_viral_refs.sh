#!/bin/bash
#SBATCH --job-name=dl_genomes
#SBATCH --output=logs/dl_genomes_%A.out
#SBATCH --error=logs/dl_genomes_%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

# -------- CONFIG --------
WORKDIR=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets
SPECIES_TAXID_FILE=$WORKDIR/viral_training_set/unique_species_taxids.tsv
OUTDIR=$WORKDIR/ref_genomes
LOG_MISSING=$OUTDIR/missing_genomes.log
GENOME_LIST=$OUTDIR/genome_list_training.tsv
mkdir -p "$OUTDIR"
> "$GENOME_LIST"
> "$LOG_MISSING"
# ------------------------

# Read lines into array
mapfile -t species_taxid_lines < "$SPECIES_TAXID_FILE"

for line in "${species_taxid_lines[@]}"; do
    # skip empty lines
    [[ -z "$line" ]] && continue

    # split line into species and taxid
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
#    wget -q -O "$genome_path" "$url"

    # Decompress the genome
    gunzip -f "$genome_path"
    genome_path="${genome_path%.gz}"

    # Add to genome list
    echo -e "${species}\t${genome_path}" >> "$GENOME_LIST"
done

echo "Finished downloading genomes and building genome list."
echo "Genome list saved at: $GENOME_LIST"
