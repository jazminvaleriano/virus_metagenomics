import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set working directory
os.chdir("/Users/jazminvaleriano/Library/Mobile Documents/com~apple~CloudDocs/03 UNIFR MS/00. SP25/00.MASTER_THESIS/FINAL_CHAPTERS/02_Assembles_evaluation/04_OM_SZ_EFFECT")

# Global font and style settings
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
})

# Load main data
df = pd.read_csv("results/03_assemblies_evaluation/combined_gen_recov_report.tsv", sep="\t")
genome_map = pd.read_csv("mock_ref_genomes/genome_list_full.tsv", sep="\t", header=None, names=["Species", "Path"])
genome_map["Assembly"] = genome_map["Path"].str.extract(r'([^/]+)\.fna')
assembly_to_species = dict(zip(genome_map["Assembly"], genome_map["Species"]))
df["Species"] = df["Assembly"].map(assembly_to_species)

# Reshape dataframe
df_long = df.melt(id_vars=["Assembly", "Assembler", "Species"], var_name="Dataset", value_name="Genome fraction (%)")
df_long["Genome fraction (%)"] = pd.to_numeric(df_long["Genome fraction (%)"], errors="coerce")
df_long = df_long.dropna(subset=["Genome fraction (%)", "Species"])

# Labels and filtering
label_map = {
    "full_sim_om_200000": "NanoSim CM 200K",
    "full_sim_om_100000": "NanoSim CM 100K",
    "full_sim_om_50000": "NanoSim CM 50K",
    "full_sim_om_20000": "NanoSim CM 20K",
    "full_sim_om_10000": "NanoSim CM 10K",
}
df_long = df_long[df_long["Dataset"].isin(label_map.keys())]
df_long["Label"] = df_long["Dataset"].map(label_map)
dataset_order = list(label_map.values())

# Load proportions
proportions_df = pd.read_csv("mock_ref_genomes/abundances_simulation.tsv", sep="\t", skiprows=1, names=["Species", "Proportion"])
proportions_df = proportions_df.set_index("Species")

virus_species = [
    "Swine Influenza A Virus (IAV-S)",
    "Betaarterivirus suid 1 (PRRSV-1)",
    "Betacoronavirus 1",
    "Porcine parvovirus",
    "Suid betaherpesvirus 2"
]
proportions_df["Group"] = proportions_df.index.map(lambda x: "Virus" if x in virus_species else "Bacteria")
proportions_df = proportions_df.sort_values(by=["Group", "Proportion"], ascending=[True, False])

# Plot one heatmap per assembler
for tool in df_long["Assembler"].unique():
    sub_df = df_long[df_long["Assembler"] == tool]
    heatmap_df = sub_df.pivot(index="Species", columns="Label", values="Genome fraction (%)")
    heatmap_df = heatmap_df.reindex(columns=dataset_order)

    # Reindex rows and extract proportions
    heatmap_df = heatmap_df.reindex(proportions_df.index)

    fig, ax = plt.subplots(figsize=(14, 10))

    sns.heatmap(
        heatmap_df,
        annot=True,
        fmt=".2f",
        cmap="YlGnBu",
        cbar_kws={"label": "Genome fraction (%)"},
        linewidths=0.5,
        linecolor="white",
        ax=ax
    )

    # Add separation line between bacteria and viruses
    bacteria_count = len(proportions_df[proportions_df["Group"] == "Bacteria"])
    ax.hlines(bacteria_count, *ax.get_xlim(), colors="black", linestyles="--", linewidth=1)

    # Add proportion values as annotations on the left
    y_labels = heatmap_df.index.tolist()
    proportions = proportions_df.loc[heatmap_df.index]["Proportion"]

    ax.set_yticklabels([
        f"{species} ({proportions[species]:.1f}%)" for species in y_labels
    ], rotation=0)

    ax.set_title(f"Genome Recovery per Species — {tool}", pad=20)
    ax.set_xlabel("")
    ax.set_ylabel("Species (Proportion in Sample)")

    plt.tight_layout()
    plt.savefig(f"figures/Genomes_rec_with_proportions_{tool}.pdf", bbox_inches="tight")
    plt.savefig(f"figures/Genomes_rec_with_proportions_{tool}.svg", bbox_inches="tight")
    plt.show()
