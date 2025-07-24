#!/bin/bash

chmod +x scripts/*

# Create logs directory
mkdir -p logs

# 1. Step 01
jid1=$(sbatch --job-name=step01 --output=logs/step01.out scripts/01_kraken_viral.sh | awk '{print $4}')

# 2. Step 02
jid2=$(sbatch --job-name=step02 --output=logs/step02.out scripts/02_centrifuger_viral.sh | awk '{print $4}')

# 3. Step 03
jid3=$(sbatch --dependency=afterok:$jid1:$jid2 --job-name=step03 --output=logs/step03.out scripts/03_get_unclassified_both.sh | awk '{print $4}')

# 4. Step 04
jid4=$(sbatch --dependency=afterok:$jid3 --job-name=step04 --output=logs/step04.out scripts/04_blast_unclassified.sh | awk '{print $4}')

# 5. Step 05
jid5=$(sbatch --dependency=afterok:$jid4 --job-name=step05 --output=logs/step05.out scripts/05_annotate_names.sh | awk '{print $4}')

# 6. Step 06
jid6=$(sbatch --dependency=afterok:$jid5 --job-name=step06 --output=logs/step06.out scripts/06_list_hits.sh | awk '{print $4}')

# 7. Step 07
jid7=$(sbatch --dependency=afterok:$jid6 --job-name=step07 --output=logs/step07.out scripts/07_prefilter_accession_map.sh | awk '{print $4}')

# 8. Step 08
jid8=$(sbatch --dependency=afterok:$jid7 --job-name=step08 --output=logs/step08.out scripts/08_Extract_ref_genomes.sh | awk '{print $4}')

# 9. Step 09
jid9=$(sbatch --dependency=afterok:$jid8 --job-name=step09 --output=logs/step09.out scripts/09_remap_viral_reads.sh | awk '{print $4}')

# 10. Step 10
jid10=$(sbatch --dependency=afterok:$jid9 --job-name=step10 --output=logs/step10.out scripts/10_mapping_metrics.sh | awk '{print $4}')

# 11. Step 11
jid11=$(sbatch --dependency=afterok:$jid10 --job-name=step11 --output=logs/step11.out scripts/11_get_true_positives.sh | awk '{print $4}')

# 12. Step 12
jid12=$(sbatch --dependency=afterok:$jid11 --job-name=step12 --output=logs/step12.out scripts/12_annotate_taxNames.sh | awk '{print $4}')

# 12. Step 13
jid13=$(sbatch --dependency=afterok:$jid12 --job-name=step12 --output=logs/step12.out scripts/13_plot_results.sh | awk '{print $4}')

echo "Pipeline submitted! Final job ID: $jid13"
