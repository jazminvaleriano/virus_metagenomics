library(vegan)
library(dplyr)
library(tibble)
library(ggplot2)
library(gt)

# ─────────────────────────────────────────────
#         MARGINAL PERMANOVA
#       & PAIRWISE COMPARISONS
# ─────────────────────────────────────────────

# Read metadata
meta <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# Remove control row
meta <- meta %>% filter(Farm != "Control")

# Set sample IDs as rownames
rownames(meta) <- meta$Sample

# Keep only relevant columns
meta_subset <- meta %>%
  select(Symptoms, Farm, Location, Age_bracket,Virus)

# Ensure they are factors
meta_subset <- meta_subset %>%
  mutate(
    Symptoms = factor(Symptoms),
    Farm = factor(Farm),
    Location = factor(Location),
    Age_bracket = factor(Age_bracket),
    Virus = factor(Virus),
    Virus_Symptoms = interaction(Virus, Symptoms, drop = TRUE)
  )

# Load beta diversity distance matrix
beta_dist <- as.matrix(read.table("02_diversity_analysis/beta_diversity_matrix.txt", header = T, row.names = 1, sep='\t', check.names = FALSE))

# Remove the control
beta_dist <- beta_dist[!rownames(beta_dist) %in% "barcode22", !colnames(beta_dist) %in% "barcode22"]

# Convert to distance object
beta_dist <- as.dist(as.matrix(beta_dist))

# Match order
sample_names <- labels(beta_dist)

meta_subset <- meta_subset %>%
  filter(rownames(meta_subset) %in% sample_names) %>%
  .[sample_names, ]

# ─────────────────────────────────────────────
#         MARGINAL PERMANOVA
# ─────────────────────────────────────────────

current_vars <- c("Farm", "Virus_Symptoms")


# Create the formula
full_formula <- as.formula(paste("beta_dist ~", paste(current_vars, collapse = " + ")))

# Run marginal PERMANOVA (type III)
marginal_model <- adonis2(
  full_formula,
  data = meta_subset,
  permutations = 999,
  method = "bray",
  by = "margin"  # ← key to marginal testing
)

# View results
print(marginal_model)

# Convert result to clean data frame
marginal_df <- as.data.frame(marginal_model) %>%
  rownames_to_column("Variable") %>%
  mutate(Significance = cut(`Pr(>F)`,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                            labels = c("***", "**", "*", ".", "ns")))

# Format into pretty table with gt
marginal_df %>%
  gt() %>%
  tab_header(title = md("**Marginal PERMANOVA**")) %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  fmt_integer(columns = "Df") %>%
  cols_align(align = "center", columns = where(is.numeric) )%>%
  cols_label(
    `SumOfSqs` = "Sum of Squares",
    `R2` = "R²",
    `F` = "F-statistic",
    `Pr(>F)` = "p-value"
  ) %>%
  tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns ≥0.1",
    locations = cells_column_labels(columns = Significance)
  )


#-------------------------------------------------
# Virus PAIRWISE
#-------------------------------------------------

# 1. Count group sizes
group_counts <- table(meta_subset$Virus)

# 2. Keep only groups with at least 2 samples
valid_groups <- names(group_counts[group_counts >= 2])

# 3. Create all pairwise combinations
pairwise_combos <- combn(valid_groups, 2, simplify = FALSE)

# 4. Run PERMANOVA for each pair
pairwise_results <- lapply(pairwise_combos, function(groups) {
  # Subset metadata
  sub_meta <- meta_subset %>% filter(Virus %in% groups)
  sample_ids <- rownames(sub_meta)
  
  # Subset distance matrix
  sub_dist <- as.dist(as.matrix(beta_dist)[sample_ids, sample_ids])
  
  # Run PERMANOVA
  model <- adonis2(sub_dist ~ Virus, data = sub_meta, permutations = 999, method = "bray")
  
  # Return stats
  data.frame(
    Group1 = groups[1],
    Group2 = groups[2],
    R2 = model$R2[1],
    F = model$F[1],
    p_value = model$`Pr(>F)`[1]
  )
})

# 5. Combine and annotate significance
pairwise_virus_df <- do.call(rbind, pairwise_results)

pairwise_virus_df <- pairwise_virus_df %>%
  mutate(Significance = cut(p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                            labels = c("***", "**", "*", ".", "ns")))

# 6. View 
print(pairwise_virus_df)

pairwise_virus_df %>%
  gt() %>%
  tab_header(title = md("**Pairwise PERMANOVA**")) %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  cols_align(align = "center", columns = where(is.numeric) )%>%
  cols_label(
    `R2` = "R²"
  ) %>%
  tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns ≥0.1",
    locations = cells_column_labels(columns = Significance)
  )
#-------------------------------------------------
# Virus_Symptoms PAIRWISE
#-------------------------------------------------

# 1. Count group sizes
group_counts <- table(meta_subset$Virus_Symptoms)

# 2. Keep only groups with at least 2 samples
valid_groups <- names(group_counts[group_counts >= 2])

# 3. Create all pairwise combinations
pairwise_combos <- combn(valid_groups, 2, simplify = FALSE)

# 4. Run PERMANOVA for each pair
pairwise_results <- lapply(pairwise_combos, function(groups) {
  # Subset metadata
  sub_meta <- meta_subset %>% filter(Virus_Symptoms %in% groups)
  sample_ids <- rownames(sub_meta)
  
  # Subset distance matrix
  sub_dist <- as.dist(as.matrix(beta_dist)[sample_ids, sample_ids])
  
  # Run PERMANOVA
  model <- adonis2(sub_dist ~ Virus_Symptoms, data = sub_meta, permutations = 999, method = "bray")
  
  # Return stats
  data.frame(
    Group1 = groups[1],
    Group2 = groups[2],
    R2 = model$R2[1],
    F = model$F[1],
    p_value = model$`Pr(>F)`[1]
  )
})

# 5. Combine and annotate significance
pairwise_df <- do.call(rbind, pairwise_results)

pairwise_df <- pairwise_df %>%
  mutate(Significance = cut(p_value,
                            breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                            labels = c("***", "**", "*", ".", "ns")))


pairwise_df %>%
  gt() %>%
  tab_header(title = md("**Pairwise PERMANOVA**")) %>%
  fmt_number(columns = where(is.numeric), decimals = 3) %>%
  cols_align(align = "center", columns = where(is.numeric) )%>%
  cols_label(
    `R2` = "R²"
  ) %>%
  tab_footnote(
    footnote = "Significance codes: *** <0.001, ** <0.01, * <0.05, . <0.1, ns ≥0.1",
    locations = cells_column_labels(columns = Significance)
  )

