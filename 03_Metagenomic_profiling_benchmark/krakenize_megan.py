import pandas as pd
import re
import sys

##CONFIG=============================
# Path to input file
input_file = "results/06_diamond-megan/final_node_with_path.tsv"
# Select one sample column to create report
sample_column = sys.argv[1] 
#====================================

# Load the merged table
df = pd.read_csv(input_file, sep="\t")


# Load the merged table
df = pd.read_csv(input_file, sep="\t")

# Function to extract the last rank code from the path
def get_last_rank(tax_path):
    if pd.isna(tax_path) or tax_path == "path":
        return 'U'
    ranks = re.findall(r'\[([A-Z])\]', str(tax_path))
    if ranks:
        return ranks[-1]
    else:
        return 'U'

# Annotate df with tax_rank column
df['tax_rank'] = df['path'].apply(get_last_rank)

# Filter rows with non-zero counts in the selected sample
filtered_df = df[df[sample_column].astype(float) > 0].copy()


# Keep only relevant columns
filtered_df = filtered_df[["taxid", "path", "tax_rank", sample_column]]

# Create accumulated_reads column
filtered_df['accumulated_reads'] = filtered_df[sample_column]

# Define rank hierarchy in reverse order
rank_order = ['K', 'D', 'P', 'C', 'O', 'F', 'G', 'S']

# Sort respecting rank hierarchy reverse
filtered_df['rank_order'] = filtered_df['tax_rank'].apply(lambda x: rank_order.index(x) if x in rank_order else 99)
filtered_df = filtered_df.sort_values(by='rank_order', ascending=True).reset_index(drop=True)

# Update accumulated_reads in reverse order
for rank in rank_order:
    for idx, row in filtered_df.iterrows():
        if row['tax_rank'] == rank:
            current_path = str(row['path']) if pd.notna(row['path']) else ''
            lower_rows = filtered_df[
                (filtered_df['path'].str.startswith(current_path, na=False)) &
                (filtered_df['path'] != current_path)
            ]
            accumulated_sum = lower_rows['accumulated_reads'].sum() + row[sample_column]
            filtered_df.at[idx, 'accumulated_reads'] = accumulated_sum


# Add percentage column
total_reads = filtered_df[sample_column].sum()
filtered_df['percentage'] = 100 * filtered_df['accumulated_reads'] / total_reads

# Drop helper column
filtered_df = filtered_df.drop(columns=['rank_order'])

# Save annotated and filtered table
output_file = input_file.replace("final_node_with_path.tsv", f"{sample_column}_filtered_with_rank.tsv")
filtered_df.to_csv(output_file, sep="\t", index=False)

print(f"Saved filtered and annotated file to {output_file}")

total_reads = filtered_df[sample_column].sum()
print(total_reads)

def extract_last_taxon_name(path):
    if pd.isna(path) or not isinstance(path, str) or ";" not in path:
        return "unclassified"
    last_element = path.split(";")[-2].strip() if path.strip().endswith(";") else path.split(";")[-1].strip()
    name = re.sub(r"\[.\]\s*", "", last_element)
    return name if name else "unclassified"

filtered_df["indented_name"] = filtered_df["path"].apply(extract_last_taxon_name)
filtered_df

report = filtered_df[["percentage", "accumulated_reads", sample_column, "tax_rank", "taxid", "indented_name"]]
output_file = input_file.replace("final_node_with_path.tsv", f"{sample_column}_kraken_style_report.tsv")
report.to_csv(output_file, sep="\t", header=False, index=False, float_format="%.2f")
print(f"Kraken-style report saved to {output_file}")

