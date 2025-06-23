#!/bin/bash
#SBATCH --job-name=centrifuge
#SBATCH --output=logs/centrifuge_dl%A.out
#SBATCH --error=logs/centrifuge_dl%A.err
#SBATCH --partition=epyc2
#SBATCH --qos=job_cpu
#SBATCH --time=2-12:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=512G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jazmin.valerianosaenz@students.unibe.ch

# This script downloads the necessary sequences to build a centrifuge db comparable to Kraken's "standard database"

module load Anaconda3
eval "$(conda shell.bash hook)"
conda activate centrifuge_env

# Centrifuge files
CENTRIFUGE_DB="databases/centrifuge_std"

#mkdir -p "$CENTRIFUGE_DB"
cd $CENTRIFUGE_DB

# Download sequences for archaea, bacterial (viral will fail with this tool, but we generate the summary): 
#centrifuge-download -o library -m -d "archaea,bacteria" refseq > seqid2taxid.map

# Adding human
#centrifuge-download -o library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' refseq >> seqid2taxid.map

# Adding viral  
SUMMARY_FILE="library/viral/assembly_summary_filtered.txt"
echo "Processing viral genomes from: $SUMMARY_FILE"

awk -F '\t' '!/^#|^$/ {print $1, $6, $20}' "$SUMMARY_FILE" | while read -r accession taxid ftp_path; do
  base=$(basename "$ftp_path")
  fna_url="${ftp_path}/${base}_genomic.fna.gz"
  fna_path="library/viral/${base}_genomic.fna.gz"

  echo "Downloading $fna_url"
  wget -q "$fna_url" -O "$fna_path"
  
  if [[ $? -eq 0 ]]; then
    echo "Downloaded: $base"
    zgrep '^>' "$fna_path" | cut -d ' ' -f1 | sed 's/^>//;s/:.*//' > seqids.txt
    awk -v taxid="$taxid" '{print $1"\t"taxid}' seqids.txt >> seqid2taxid.map
    rm seqids.txt
  else
    echo "Failed to download: $fna_url" >&2
  fi
done
conda deactivate
