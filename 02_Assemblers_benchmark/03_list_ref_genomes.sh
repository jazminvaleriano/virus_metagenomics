#!/bin/bash

ABUNDANCES_SIM_FILE=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/full_training_set/abundances_simulation.tsv
GENOMES_MAP=/storage/research/vetsuisse_ivi/jvaleriano/01_Simulated_datasets/ref_genomes/genome_list_training_full.tsv

OUT_FILE=sim_DS_genomes_list.txt
OUT_DIR=simulated_DS_references

mkdir -p $OUT_DIR

# Empty the output file
> "$OUT_FILE"

# Skip header and get species names
tail -n +2 "$ABUNDANCES_SIM_FILE" | cut -f1 | while read -r species; do
    # Find the corresponding genome path
    path=$(awk -v sp="$species" -F'\t' '$1 == sp {print $2}' "$GENOMES_MAP")

    if [[ -n "$path" ]]; then
        echo "$path" >> "$OUT_FILE"
    else
        echo "Warning: No path found for $species"
    fi
done

# Create symlinks for each path
while read -r filepath; do
    ln -s "$filepath" "$OUT_DIR"
done < "$OUT_FILE"