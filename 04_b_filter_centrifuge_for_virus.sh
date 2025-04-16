#!/bin/bash

mkdir -p results/centrifuge/viruses_only

for file in results/centrifuge/*.centrifuge.report; do
    out="results/centrifuge/viruses_only/$(basename "$file" .centrifuge.report).viruses_only.tsv"
    awk 'NR==1 || tolower($1) ~ /virus/' "$file" > "$out"
done
