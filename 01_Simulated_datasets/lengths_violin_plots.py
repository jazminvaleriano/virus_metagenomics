import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# File paths
files = {
    "sw-eSISPA-ONT": "read_lengths/real_DS_filt.tsv",
    "ZM Complete": "read_lengths/Zymo.tsv",
    "ZM Short 200-20k bp": "read_lengths/Zymo_SHORT.tsv",
    "ZM Long 2k-50k bp": "read_lengths/Zymo_LONG.tsv",
    "ZM Short adjusted": "read_lengths/Zymo_SHORT.tsv",
    "ZM Long adjusted": "read_lengths/Zymo_LONG.tsv",
    "Full simulated": "read_lengths/Full_sim.tsv",
    "Full simulated_pretrained": "read_lengths/Full_sim_pretrained.tsv"
}

# Define barcode filters for SHORT/LONG variants
barcode_filter = {
    "ZM Short 200-20k bp": "02_Zymo-200K_SHORT",
    "ZM Long 2k-50k bp": "02_Zymo-200K_LONG",
    "ZM Short adjusted": "03_Zymo-200K_length_adjusted_SHORT",
    "ZM Long adjusted": "03_Zymo-200K_length_adjusted_LONG",
    "Full simulated": "full_simulated_dataset_200000",
    "Full simulated_pretrained": "full_simulated_dataset_pt_200000"
}

# Initialize list for storing all data
all_samples = []

# Load real_DS, Zymo Original and with subsampling
real_ds = pd.read_csv(files["sw-eSISPA-ONT"], sep="\t")
real_ds_sample = real_ds.sample(n=200_000, random_state=42)
real_ds_sample["sample"] = "sw-eSISPA-ONT"
all_samples.append(real_ds_sample[["sample", "read_length","log_length"]])

zymo_orig = pd.read_csv(files["ZM Complete"], sep="\t")
zymo_sample = zymo_orig.sample(n=200_000, random_state=42)
zymo_sample["sample"] = "ZM Complete"
all_samples.append(zymo_sample[["sample", "read_length","log_length"]])

# Load and filter the sample-specific samples
for sample_name in ["ZM Short 200-20k bp", "ZM Long 2k-50k bp", "ZM Short adjusted", "ZM Long adjusted","Full simulated","Full simulated_pretrained"]:
    df = pd.read_csv(files[sample_name], sep="\t")
    filtered = df[df["barcode"] == barcode_filter[sample_name]].copy()
    filtered["sample"] = sample_name
    all_samples.append(filtered[["sample","read_length", "log_length"]])

# Combine into one dataframe
final_df = pd.concat(all_samples, ignore_index=True)

# Preview
print(final_df["sample"].value_counts())
print(final_df.head())


# Set global style
sns.set(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Define sample order
sample_order = [
    "sw-eSISPA-ONT",
    "ZM Complete",
    "ZM Short 200-20k bp",
    "ZM Short adjusted",
    "ZM Long 2k-50k bp",
    "ZM Long adjusted",
    "Full simulated",
    "Full simulated_pretrained"
]

# Define color map: use identical colors for paired samples
color_map = {
    "sw-eSISPA-ONT": "#1b9e77",          # teal green
    "ZM Complete": "#7570b3",    # purple
    "ZM Short 200-20k bp": "#d95f02",            # orange
    "ZM Short adjusted": "#d95f02", # orange
    "ZM Long 2k-50k bp": "#e7298a",             # pink
    "ZM Long adjusted": "#e7298a",   # pink
    "Full simulated": "#438cd1",
    "Full simulated_pretrained": "#438cd1"

}
# Create the plots
plt.figure(figsize=(14, 6))
sns.violinplot(
    data=final_df,
    x="sample",
    y="read_length",
    hue="sample",
    palette=color_map,
    order=sample_order,
    inner="quartile",
    bw=0.4,
    cut=0,
    linewidth=1,
    legend=False
)

# Axis formatting
plt.ylim(0, 20000)
plt.ylabel("log Length (bp)", labelpad=10)
plt.xlabel("")
plt.xticks(rotation=0)
plt.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()

# Save figures
# plt.savefig("log_lengths_violin_pub.pdf", bbox_inches="tight")
# plt.savefig("log_lengths_violin_pub.svg", bbox_inches="tight")
# plt.savefig("log_lengths_violin_pub.png", dpi=600)

plt.show()

# Create the plots
plt.figure(figsize=(14, 6))
sns.violinplot(
    data=final_df,
    x="sample",
    y="log_length",
    hue="sample",
    palette=color_map,
    order=sample_order,
    inner="quartile",
    bw=0.4,
    cut=0,
    linewidth=1,
    legend=False
)

# Axis formatting
plt.ylim(0, 6)
plt.ylabel("log Length (bp)", labelpad=10)
plt.xlabel("")
plt.xticks(rotation=0)
plt.grid(True, axis='y', linestyle='--', linewidth=0.5, alpha=0.7)
plt.tight_layout()

# Save figures
# plt.savefig("log_lengths_violin_pub.pdf", bbox_inches="tight")
# plt.savefig("log_lengths_violin_pub.svg", bbox_inches="tight")
# plt.savefig("log_lengths_violin_pub.png", dpi=600)

plt.show()