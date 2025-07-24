#!/bin/bash
#SBATCH --job-name=att_names
#SBATCH --error=logs/att_names_%A.err

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

INPUT_DIR=04_positive_results
TAXONKIT_OUTDIR=04_positive_results
FASTA_DB=/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/refseq_viral/viral.1.1.genomic.fna
ACCESSION_NAMES="$TAXONKIT_OUTDIR/accession_names.tsv"

mkdir -p "$TAXONKIT_OUTDIR"

# Extract accession names from FASTA only once
if [[ ! -f "$ACCESSION_NAMES" ]]; then
    echo "Extracting accession names from FASTA..."
    grep "^>" "$FASTA_DB" | sed 's/^>//' | awk '{acc=$1; $1=""; print acc "\t" substr($0,2)}' > "$ACCESSION_NAMES"
fi

# Loop over metric files
for METRICS_FILE in "$INPUT_DIR"/*_viral_metrics.tsv; do
    SAMPLE=$(basename "$METRICS_FILE" _viral_metrics.tsv)
    TAXIDS="$TAXONKIT_OUTDIR/${SAMPLE}_taxids.tsv"
    LINEAGE="$TAXONKIT_OUTDIR/${SAMPLE}_taxids_lineage.tsv"
    FINAL="$TAXONKIT_OUTDIR/${SAMPLE}_viral_metrics_annotated.tsv"

    echo "Annotating $SAMPLE..."

    # Extract taxIDs from metrics
    cut -f2 "$METRICS_FILE" | tail -n +2 | sort -u > "$TAXIDS"

    # Get lineage from taxonkit
    taxonkit lineage "$TAXIDS" > "$LINEAGE"

    # Merge in Python
    python3 - <<EOF
import pandas as pd

metrics = pd.read_csv("$METRICS_FILE", sep="\t")
lineage = pd.read_csv("$LINEAGE", sep="\t", names=["Taxon", "Lineage"], dtype=str)
names = pd.read_csv("$ACCESSION_NAMES", sep="\t", names=["Accession", "AccessionName"], dtype=str)

# Extract species name from lineage
lineage["TaxonName"] = lineage["Lineage"].str.split(";").str[-1].str.strip()

# Clean merge keys
metrics["Taxon"] = metrics["Taxon"].astype(str).str.strip()
lineage["Taxon"] = lineage["Taxon"].astype(str).str.strip()
metrics["Accession"] = metrics["Accession"].astype(str).str.strip()
names["Accession"] = names["Accession"].astype(str).str.strip()

# Merge all
merged = metrics.merge(lineage, how="left", on="Taxon")
merged = merged.merge(names, how="left", on="Accession")

# Save output
merged.to_csv("$FINAL", sep="\t", index=False)
EOF

    echo "Annotated: $FINAL"
    mv $FINAL $METRICS_FILE

done

# Cleanup
rm ${TAXONKIT_OUTDIR}/*taxids*
rm ${TAXONKIT_OUTDIR}/accession_names.tsv


echo "All annotation steps completed."
