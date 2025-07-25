#!/bin/bash 

CONCAT_FILE=results/species_concatenated_precision.csv

echo -e "Classifier,True Positives,False Positives,False Negatives,Precision,Recall,F1,F0.5,threshold,dataset" > $CONCAT_FILE

# Define only the folders you want
for folder in "100K_ZM Long (length-adjusted)" "20K_ZM Long (length-adjusted)" "5K_ZM Long (length-adjusted)"; do
    for f in "$folder"/*-Species-Precision-Recall*; do
        dataset=$(basename "$f")
        dataset="${dataset%%-Species-Precision-Recall*}"
        echo "Processing $dataset"

        tail -n +2 "$f" | awk -v ds="$dataset" -F',' 'BEGIN{OFS=","} {print $0, ds}' >> $CONCAT_FILE
    done
done
