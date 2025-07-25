import os
os.chdir("/Users/jazminvaleriano/Library/Mobile Documents/com~apple~CloudDocs/03 UNIFR MS/00. SP25/00.MASTER_THESIS/FINAL_CHAPTERS/03_Bacterial_mNGS_validation_LONG")
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv("results/concatenated_precision.csv")

dataset_order = [
#    "200K_ZM Long (length-adjusted)",
    "5K_ZM Long (length-adjusted)",
    "20K_ZM Long (length-adjusted)",
    "100K_ZM Long (length-adjusted)",
#    "10K_ZM Long (length-adjusted)",
]
custom_labels = [#"200K ZM",
                "5K ZM", 
                 "20K ZM",
                 "100K ZM",
                 #"10K ZM",
         ]
df['dataset'] = pd.Categorical(df['dataset'], categories=dataset_order, ordered=True)

pal_class = sns.color_palette(['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02'])

# Melt the dataframe to have Metric as a variable
df_melted = df.melt(
    id_vars=['dataset', 'Classifier', 'threshold'],
    value_vars=['Precision', 'Recall', 'F1'],
    var_name='Metric',
    value_name='Score'
)

sns.set_theme(style="whitegrid", font_scale=1)

g = sns.catplot(
    data=df_melted,
    kind="bar",
    x="dataset",
    y="Score",
    hue="Classifier",
    col="threshold",
    row="Metric",
    palette=pal_class,
    height=4,
    aspect=1.2,
    dodge=True,
    sharey=False,
    legend_out=True
)

for ax in g.axes.flatten():
    ax.set_xticks(range(len(custom_labels)))
    ax.set_xticklabels(custom_labels) #rotation=45, ha="right")
    ax.grid(color="gray", linestyle="--", linewidth=0.5, axis="y")
    ax.set_xlabel("")


sns.despine(trim=True)

g.set_titles(row_template="{row_name}", col_template="Threshold = {col_name}")
g.savefig("results/precision_recall_f1_plot.png", dpi=300, bbox_inches='tight')
plt.show()

heat_data = df_melted[df_melted["Metric"] == "F1"].pivot_table(
    index="dataset", columns="Classifier", values="Score"
)

sns.heatmap(heat_data, annot=True, fmt=".2f", cmap="viridis", cbar_kws={'label': 'F1 Score'})
plt.title("F1 Score by Dataset and Classifier")
plt.ylabel("Dataset")
plt.xlabel("Classifier")
plt.tight_layout()

plt.show()


# Trying faceted point plot
import seaborn as sns
import matplotlib.pyplot as plt

# Update threshold to string or nicely formatted float for x-axis
df_melted["threshold_str"] = df_melted["threshold"].astype(str)

# Order dataset columns
dataset_order = [
    "5K_ZM Long (length-adjusted)",
    "20K_ZM Long (length-adjusted)",
    "100K_ZM Long (length-adjusted)"
]
custom_labels = ["5K ZM", "20K ZM", "100K ZM"]
df_melted['dataset'] = pd.Categorical(df_melted['dataset'], categories=dataset_order, ordered=True)

# Plot
g = sns.catplot(
    data=df_melted,
    kind="point",
    x="threshold_str",
    y="Score",
    hue="Classifier",
    col="dataset",
    row="Metric",
    palette=pal_class,
    height=4,
    aspect=1.2,
    dodge=0.3,
    markers="o",
    linestyles="-"
)

# Set axis and column titles
for ax, label in zip(g.axes[0], custom_labels):
    ax.set_title(label)

for ax in g.axes.flatten():
    ax.grid(color="gray", linestyle="--", linewidth=0.5, axis="y")
    ax.set_xlabel("Threshold")

sns.despine(trim=True)
g.set_axis_labels("Threshold", "Score")
g.set_titles(row_template="{row_name}", col_template="{col_name}")
g.savefig("results/precision_recall_f1_pointplot_by_size.svg", dpi=300, bbox_inches="tight")

g.set_titles(row_template="{row_name}", col_template="Threshold = {col_name}")
g.savefig("results/precision_recall_f1_point.png", dpi=300, bbox_inches='tight')

plt.tight_layout()
plt.show()
