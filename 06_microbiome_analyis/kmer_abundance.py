import sys
import pandas as pd

# ---------------- Parse Arguments ----------------
file_path = sys.argv[1]
taxonomic_level = sys.argv[2]
min_threshold = int(sys.argv[3])
output_file = f"{file_path}_abundance.tsv"

# ---------------- Load and Process ----------------
df = pd.read_table(file_path, header=None, sep='\t')
lineage_col = df.columns[-1]  # I added the lineage annotation as a last step of Kraken classification

# Filter by taxonomic rank
df = df[df[5] == taxonomic_level]

# Filter by distinct minimizer threshold (leave out taxa that were classified based on a very short sequence (even if many reads matched it))
df = df[df[4] >= min_threshold]

# Keep only Bacteria and Archaea, filter host and human contamination out
df = df[
    (df[lineage_col].str.contains("Bacteria|Archaea", case=False, na=False)) &
    (~df[lineage_col].str.contains("Viruses|Eukaryota|Homo sapiens|Sus scrofa", case=False, na=False))
]

# Calculate relative abundance based on mapped minimizers
df['rel_abundance'] = df[3] / df[3].sum() 

# Save result
df.to_csv(output_file, sep='\t', index=False, header=False)

print(f"Output saved to: {output_file}")
