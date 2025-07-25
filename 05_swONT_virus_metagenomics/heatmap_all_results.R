library(tidyverse)
library(pheatmap)
library(viridis)

# 1. Read and clean data
df <- read_tsv("04_positive_results/combined_report.tsv", 
               col_select = c("Sample","Lineage", "R", "Farm", "Symptoms", "Pass")) %>%
  filter(Pass == TRUE)

# Add column for genus
df <- df %>%
  mutate(Genus = str_split(Lineage, ";", simplify = TRUE)[, 9])


# 2. Filter unclassified or unwanted genera
excluded_genera <- c(
  "Gemykibivirus anima1",
  "Gemykrogvirus carib1",
  "Porcine stool-associated circular virus/BEL/15V010",
  "Fur seal faeces associated circular DNA virus"
)
df <- df %>% filter(!Genus %in% excluded_genera)


# 4. Pivot to wide format for heatmap
heatmap_data <- df %>%
  select(Sample, Genus, R) %>%
  group_by(Sample, Genus) %>%
  summarise(R = max(R), .groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = R, values_fill = 0)

# 5. Convert to numeric matrix
mat <- as.data.frame(heatmap_data)
rownames(mat) <- mat$Sample
mat <- mat %>% select(-Sample)
mat <- as.matrix(mat)

# 6. Metadata annotations for rows
annotation_row <- df %>%
  distinct(Sample, Farm, Symptoms) %>%
  column_to_rownames("Sample") %>%
  mutate(
    Farm = factor(Farm),
    Symptoms = factor(Symptoms)
  )

# 7. Define annotation colors
annotation_colors <- list(
  Symptoms = c(
    "With Symptoms" = "#e41a1c",
    "No Symptoms" = "#4daf4a"
  ),
  Farm = c(
    "Farm 3" = "#1b9e77",
    "Farm 4" = "#d95f02",
    "Farm 5" = "#7570b3",
    "Farm 9" = "#e7298a",
    "Farm 10" = "#66a61e",
    "Farm 18" = "#e6ab02",
    "Farm 33" = "#a6761d",
    "Farm 74" = "#666666"
  )
)

# 8. Create a matrix of R values only if > 0.2
numbers_matrix <- ifelse(mat > 0.2, sprintf("%.1f", mat), "")

# 9. Plot the heatmap
pheatmap(
  mat,
  annotation_row = annotation_row,
  annotation_colors = annotation_colors,
#  display_numbers = numbers_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "R values per Sample and Genus",
  color = colorRampPalette(c("white", "#FDBB84", "#E34A33", "#B30000"))(100),
  fontsize = 11,
  fontsize_row = 9,
  fontsize_col = 10,
  border_color = "grey90",
  angle_col = 90
)
