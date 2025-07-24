#!/usr/bin/env python3
import sys
import pandas as pd
from collections import defaultdict

# Get positional arguments
if len(sys.argv) != 3:
    print("Usage: python krakenize_diamond.py <input_file> <output_file>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

# Load TSV file
df = pd.read_csv(input_file, sep="\t", header=None,
                 names=["read_id", "taxon_id", "evalue", "path"])

df["path"] = df["path"].fillna("")
summary = df.groupby("path").size().reset_index(name="reads")
summary_taxids = df.drop_duplicates(subset="path")[["path", "taxon_id"]]
summary = summary.merge(summary_taxids, on="path", how="left")

def get_rank(path):
    if path == "": return "U"
    levels = [p.strip() for p in path.split(";") if p.strip()]
    return ["R", "R1", "D", "P", "C", "O", "F", "G", "S", "S1"][min(len(levels)-1, 9)]

summary["tax_rank"] = summary["path"].apply(get_rank)

TOTAL_READS = len(df)
path_to_reads = dict(zip(summary["path"], summary["reads"]))
sorted_paths = sorted(path_to_reads.keys(), key=lambda x: x.count(";"))
cumulative_reads = defaultdict(int)

for path in sorted_paths:
    cumulative_reads[path] = path_to_reads[path]

for parent in sorted_paths:
    for child in sorted_paths:
        if child != parent and child.startswith(parent + ";"):
            cumulative_reads[parent] += path_to_reads[child]

summary["cumulative_reads"] = summary["path"].map(cumulative_reads)
summary["percentage"] = 100 * summary["cumulative_reads"] / TOTAL_READS

summary["indented_path"] = summary.apply(lambda row: "  " * row["path"].count(";") + row["path"].split(";")[-1] if row["path"] else "unclassified", axis=1)
report = summary[["percentage", "cumulative_reads", "reads", "tax_rank", "taxon_id", "indented_path"]]
report.to_csv(output_file, sep="\t", header=False, index=False, float_format="%.2f")
print(f"Report saved to: {output_file}")