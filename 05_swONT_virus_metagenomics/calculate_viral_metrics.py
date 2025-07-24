import sys
import math
import pandas as pd

# --- Inputs ---
depth_file = sys.argv[1]               # e.g. sample_depth.txt
accession_stats_file = sys.argv[2]     # e.g. sample_accession_stats.tsv
taxon_map_file = sys.argv[3]           # e.g. sample_accession_taxon_map.tsv
out_file = sys.argv[4]                 # output TSV

# --- Load data ---
# Depth info
depth = pd.read_csv(depth_file, sep="\t", names=["Accession", "Position", "Depth"])
depth["Covered"] = depth["Depth"] >= 1
covered_bases = depth.groupby("Accession")["Covered"].sum().reset_index(name="CoveredBases")

# Stats
stats = pd.read_csv(accession_stats_file, sep="\t")

# ❗ Remove accessions with no mapped reads
stats = stats[stats["MappedReads"] > 0].copy()

# Taxon info
taxa = pd.read_csv(taxon_map_file, sep="\t", header=None, usecols=[1,2], names=["Accession", "Taxon"])

# Merge all
df = stats.merge(covered_bases, on="Accession", how="left")
df["CoveredBases"] = df["CoveredBases"].fillna(0)
df = df.merge(taxa, on="Accession", how="left")

# --- Calculate metrics ---
df["Co"] = df["CoveredBases"] / df["GenomeLength"]
df["D"] = (df["MappedReads"] * df["AvgReadLength"]) / df["GenomeLength"]
df["Ce"] = 1 - df["D"].apply(lambda d: math.exp(-d) if d < 700 else 0)  # avoid overflow
df["R"] = df["Co"] / df["Ce"]
df["Pass"] = (df["R"] >= 0.3) & (df["Co"] >= 0.1)

# --- Output ---
df[["Accession", "Taxon", "GenomeLength", "MappedReads", "AvgReadLength", "CoveredBases", "Co", "Ce", "R", "Pass"]] \
  .to_csv(out_file, sep="\t", index=False)

print(f"✅ Metrics with taxon names written to {out_file}")
