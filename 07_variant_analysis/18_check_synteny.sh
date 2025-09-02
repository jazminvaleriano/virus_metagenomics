#!/bin/bash
#SBATCH --job-name=synteny
#SBATCH --output=logs/synteny_%A_%a.out
#SBATCH --error=logs/synteny_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# Load environment
module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate metagenomics_env

GENOMES=(/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references/parvovirus1/parvovirus1.fasta \
/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references/Porcine_Bocavirus/Porcine_Bocavirus.fasta \
/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references/Porcine_Hokovirus/Porcine_Hokovirus.fasta \
/storage/research/vetsuisse_ivi/jvaleriano/Final_databases/references/all_parvoviridae/parvovirus3.fasta)

# Make output dir
mkdir -p 03_mummer_out_filtered
cd 03_mummer_out_filtered

# All pairwise comparisons
for ((i=0; i<${#GENOMES[@]}; i++)); do
  for ((j=i+1; j<${#GENOMES[@]}; j++)); do
    REF="${GENOMES[$i]}"
    QRY="${GENOMES[$j]}"
    PFX="$(basename "${REF%.*}")_vs_$(basename "${QRY%.*}")"

    echo ">>> Aligning $REF vs $QRY"

    # 1) Align (permissive seeding to detect any similarity)
    nucmer --maxmatch -l 5 -c 40 -p "$PFX" "$REF" "$QRY"

    # 2) Filter with mapper-like thresholds (reportable / "could recruit reads")
    #    keep alignments >=100 bp and >=80% identity; many-to-many to retain any plausible block
    delta-filter -m -i 75 -l 80 "${PFX}.delta" > "${PFX}.len100_i80.delta"
    show-coords -rclT "${PFX}.len100_i80.delta" > "${PFX}.len100_i80.coords.tsv"

    # 3) Optional stricter view for a clean figure (higher identity/length, 1-to-1)
    delta-filter -1 -i 90 -l 300 "${PFX}.delta" > "${PFX}.len300_i90.1to1.delta"
    show-coords -rclT "${PFX}.len300_i90.1to1.delta" > "${PFX}.len300_i90.1to1.coords.tsv"

    # 4) Plots
    #    Exploratory/permissive plot (may be empty if truly divergent)
    mummerplot --png --large --filter -p "${PFX}.len100_i80" "${PFX}.len100_i80.delta"
    #    Clean, stricter plot
    mummerplot --png --large --filter -p "${PFX}.len300_i90.1to1" "${PFX}.len300_i90.1to1.delta"
    mummerplot --pdf --large --filter -p "${PFX}.len300_i90.1to1" "${PFX}.len300_i90.1to1.delta"
  done
done
