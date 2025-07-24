import pandas as pd

# CONFIG
input_file = "full_training_set/species_read_counts_summary.tsv"

# Load your full species summary file
df = pd.read_csv(input_file, sep="\t", names=("Species", "Taxon", "reads", "Proportion"))

# Step 1: Stratified selection
high_abundance = df.iloc[:20]
medium_abundance = df.iloc[20:100]
low_abundance = df.iloc[100:]

# Select species per stratum
high_selected = high_abundance.head(10)
medium_selected = medium_abundance.head(5)
low_selected = low_abundance.head(5)

# Combine selected species
selected_bacteria = pd.concat([high_selected, medium_selected, low_selected]).reset_index(drop=True)

print(selected_bacteria)

# Step 3: Rescale proportions to sum 98%
total_selected_prop = selected_bacteria['Proportion'].sum()
selected_bacteria['New_Proportion'] = (selected_bacteria['Proportion'] / total_selected_prop) * 98
selected_bacteria['New_Proportion'] = selected_bacteria['New_Proportion'].round(2)

# Step 4: Add fixed viruses
virus_data = [
    ["Swine Influenza A Virus (IAV-S)", "Virus", 0.14],
    ["Betaarterivirus suid 1 (PRRSV-1)", "Virus", 0.14],
    ["Betacoronavirus 1", "Virus", 0.30],
    ["Porcine parvovirus", "Virus", 0.06],
    ["Suid betaherpesvirus 2", "Virus", 1.36]
]
virus_df = pd.DataFrame(virus_data, columns=["Species", "Taxon", "New_Proportion"])
virus_df['New_Proportion'] = virus_df['New_Proportion'].round(2)

# Step 5: Merge and export
final_df = pd.concat([
    selected_bacteria[['Species', 'Taxon', 'New_Proportion']],
    virus_df
]).reset_index(drop=True)

# Optional: Check if total is exactly 100 and print warning if not
total_prop_sum = final_df['New_Proportion'].sum()
print(f"Total proportion sum: {total_prop_sum:.2f}%")

final_df.to_csv("full_training_set/mock_community_proportions.tsv", sep="\t", index=False)

print("Mock community proportions saved as 'full_training_set/mock_community_proportions.tsv'")
