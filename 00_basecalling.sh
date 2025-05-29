#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=dorado_%j.out

POD_FILES="/storage/research/vetsuisse_ivi/jvaleriano/20250324_1831_MN34349_FAY63474_1833a9d4/pod5"
OUT_DIR="/storage/research/vetsuisse_ivi/jvaleriano/20250324_1831_MN34349_FAY63474_1833a9d4/RecalledReads_may29"

module load apptainer

apptainer run dorado.sif basecaller \
  --model ./models/dna_r10.4.1_e8.2_400bps_hac@v4.3.0 \
  --flowcell FLO-MIN114 \
  --kit SQK-NBD114-96 \
  --barcode-trim \
  --min-qscore 9 \
  $POD_FILES ./output
