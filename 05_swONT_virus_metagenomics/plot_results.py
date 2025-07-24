import sys
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Rectangle

filename=sys.argv[1]  
sample_name = os.path.basename(filename).split("_")[0]
directory= os.path.dirname(filename)

df = pd.read_csv(filename, sep="\t")

# Make sure 'Lineage' is treated as string and doesn't have trailing semicolons
df["Lineage"] = df["Lineage"].astype(str).str.strip().str.rstrip(";")

# Split lineage by semicolon and get the last two fields
df["LastTwoLineage"] = df["Lineage"].apply(
    lambda x: ";".join(x.split(";")[-2:]) if ";" in x else x
)

# Sort by R for clarity
df_sorted = df.sort_values("R", ascending=False)


##====================METRICS PLOT=================================

# Add FalsePositive flag
df["FalsePositive"] = (df["R"] < 0.3) | (df["Co"] < 0.1)
df["Label"] = df["Accession"] + " | " + df["LastTwoLineage"]

# Prepare data
heatmap_data = df.set_index("Label")[["Co", "Ce", "R"]]
mask = df.set_index("Label")["FalsePositive"].reindex(heatmap_data.index)

# Sort heatmap to group false positives
heatmap_data["FalsePositive"] = mask
heatmap_data = heatmap_data.sort_values(by="FalsePositive", ascending=False)
mask = heatmap_data.pop("FalsePositive")  # remove the column but keep the order

# Choose a better diverging colormap
cmap = sns.light_palette("green", as_cmap=True)

# Create the heatmap
fig, ax = plt.subplots(figsize=(12, max(4, len(heatmap_data) * 0.6)) )
sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".2f",
    cmap=cmap,
    linewidths=0.6,
    linecolor="white",
    cbar_kws={'label': 'Metric Value'},
    annot_kws={"fontsize": 10},
    ax=ax
)

# Add red rectangles for false positives
for y_index, label in enumerate(heatmap_data.index):
    if mask[label]:
        ax.add_patch(Rectangle((0, y_index), 3, 1, fill=False, edgecolor='red', lw=2, clip_on=False))

# Labels and title
plt.title(f"Coverage Metrics - {sample_name}", fontsize=14, fontweight="bold")
ax.set_xlabel("Metric", fontsize=13)
ax.set_ylabel("Accession ID | Taxon", fontsize=13)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=11)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, rotation=0)

# Legend / Note
plt.figtext(0.98, 0.02,
            "Red box: Potential false positive (R < 0.3 or Co < 0.1)",
            ha="right", fontsize=9, style='italic')

plt.tight_layout(pad=3.0)

plt.savefig(f"{directory}/{sample_name}_metrics_heatmap.svg", dpi=300)

##====================MAPPED READS PLOT=================================

# 1. Group by AccessionName and take the max MappedReads in case of duplicates
heatmap_reads = df.groupby("AccessionName", as_index=True)["MappedReads"].max().to_frame()

# Truncate long accession names
heatmap_reads.index = heatmap_reads.index.str.slice(0, 100) 

# Sort by MappedReads
heatmap_reads = heatmap_reads.sort_values("MappedReads", ascending=False)

# Create figure
n_rows = heatmap_reads.shape[0]
fig_height = max(3, n_rows * 0.5)  # scale dynamically depending on ds length
fig, ax = plt.subplots(figsize=(12, fig_height))

sns.heatmap(
    heatmap_reads,
    annot=True,
    fmt=".0f",
    cmap="YlOrRd",
    linewidths=0.3,
    linecolor="white",
    cbar_kws={'label': 'Mapped Reads'},
    annot_kws={"fontsize": 9},
    ax=ax
)

# Clean labels
ax.set_title(f"Mapped Reads per Viral Accession - {sample_name}", fontsize=12, fontweight="bold")
ax.set_xlabel("") 
ax.set_ylabel("Accession", fontsize=11)
ax.set_xticklabels([])
ax.set_yticklabels(ax.get_yticklabels(), fontsize=9, rotation=0)

# Fix layout
plt.tight_layout()

plt.savefig(f"{directory}/{sample_name}_mapped_reads.svg", dpi=300)
