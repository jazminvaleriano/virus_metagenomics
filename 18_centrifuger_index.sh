#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH --output=logs/centrifuger_idx%A.out
#SBATCH --error=logs/centrifuger_idx%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment (adjust if needed)
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuger_env

# -------------------- CONFIG --------------------
CENTRIFUGE_DB="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/databases/centrifuge_std"
FASTA="$CENTRIFUGE_DB/input-sequences.fna"    # Reference fasta built for CENTRIFUGE
SEQID2TAX="$CENTRIFUGE_DB/seqid2taxid.map"    # Mapping file: seqID <tab> taxID
KRAKEN_DB="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/databases/kraken2_std"
NODES="$KRAKEN_DB/taxonomy/nodes.dmp"
NAMES="$KRAKEN_DB/taxonomy/names.dmp"

CENTRIFUGER_DB="/storage/research/vetsuisse_ivi/jvaleriano/metagenomics_project/databases/centrifuger_std"        
MEMORY="400G"
# ------------------------------------------------

# Create index folder if needed
mkdir -p "$(dirname "$CENTRIFUGER_DB")"
cd $CENTRIFUGER_DB

# Build Centrifuger index
echo "Building Centrifuger index..."
centrifuger-build -r "$FASTA" \
  --conversion-table "$SEQID2TAX" \
  --taxonomy-tree "$NODES" \
  --name-table "$NAMES" \
  -t "$SLURM_CPUS_PER_TASK" \
  --build-mem "$MEMORY"

conda deactivate