import sys
from pathlib import Path

genomes_list = Path(sys.argv[1])      # e.g. SAMPLE_hit_genomes_list.tsv
fasta_file = Path(sys.argv[2])        # e.g. viral.1.1.genomic.fna
output_fasta = Path(sys.argv[3])      # e.g. SAMPLE_viral_hits.fasta

# Load desired accessions (column 2 in genomes_list)
with genomes_list.open() as f:
    wanted = set(line.strip().split("\t")[1] for line in f if line.strip())

# Search and extract matching sequences
write_seq = False
with fasta_file.open() as infile, output_fasta.open("w") as out:
    for line in infile:
        if line.startswith(">"):
            acc = line[1:].split()[0]
            write_seq = acc in wanted
        if write_seq:
            out.write(line)
