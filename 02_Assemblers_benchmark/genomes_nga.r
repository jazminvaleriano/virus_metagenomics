library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Load data
df_om <- read_tsv("04_OM_SZ_EFFECT/results/03_assemblies_evaluation/combined_nga50_report.tsv",
               col_types = cols(.default = "c"))
df_pt <- read_tsv("05_PT_SZ_EFFECT/results/03_assemblies_evaluation/combined_nga50_report.tsv",
               col_types = cols(.default = "c"))

genome_map <- read_tsv("04_OM_SZ_EFFECT/mock_ref_genomes/genome_list_full.tsv", col_names = c("Species", "Path"))
proportions_df_om <- read_tsv("04_OM_SZ_EFFECT/mock_ref_genomes/abundances_simulation.tsv", skip = 1, col_names = c("Species", "Proportion"))
lengths_df_om <- read_tsv("04_OM_SZ_EFFECT/results/03_assemblies_evaluation/genomes_length.tsv", skip = 1, col_names = c("Genome", "Length"))

# Convert assembler result columns to numeric, force invalid values to 0
assembler_cols <- grep("^full_sim_om_", names(df_om), value = TRUE)

df_om <- df_om %>%
  mutate(across(all_of(assembler_cols), ~ as.numeric(trimws(.)))) %>%
  mutate(across(all_of(assembler_cols), ~ replace_na(., 0)))

# Map Reference length
df_om <- df_om %>%
  left_join(lengths_df_om, by = c("Assembly" = "Genome")) %>%
  rename(Reference_length = Length)

# Map genome names to Species
genome_map <- genome_map %>%
  mutate(Assembly = gsub(".*/|\\.fna", "", Path))

df_om <- df_om %>%
  left_join(genome_map %>% select(Assembly, Species), by = "Assembly")

# Reshape dataframe
df_om_long <- df_om %>%
  pivot_longer(cols = all_of(assembler_cols),
               names_to = "Dataset", values_to = "NGA50")

# Labels and filtering
label_map_om <- c(
  "full_sim_om_200000" = "NanoSim CM 200K",
  "full_sim_om_100000" = "NanoSim CM 100K",
  "full_sim_om_50000" = "NanoSim CM 50K",
  "full_sim_om_20000" = "NanoSim CM 20K",
  "full_sim_om_10000" = "NanoSim CM 10K"
)

df_om_long <- df_om_long %>%
  filter(Dataset %in% names(label_map_om)) %>%
  mutate(Label = label_map_om[Dataset])

virus_species <- c(
  "Swine Influenza A Virus (IAV-S)",
  "Betaarterivirus suid 1 (PRRSV-1)",
  "Betacoronavirus 1",
  "Porcine parvovirus",
  "Suid betaherpesvirus 2"
)

df_om_long <- df_om_long %>%
  filter(Species %in% virus_species) %>%
  mutate(NGA50_Reference_Length = NGA50 / Reference_length)

df_om_filtered <- df_om_long %>%
  filter(Label %in% c("NanoSim CM 200K", "NanoSim CM 20K"))

# --- Repeat same cleaning steps for df_pt ---
assembler_cols_pt <- grep("^full_sim_pt_", names(df_pt), value = TRUE)

df_pt <- df_pt %>%
  mutate(across(all_of(assembler_cols_pt), ~ as.numeric(trimws(.)))) %>%
  mutate(across(all_of(assembler_cols_pt), ~ replace_na(., 0)))

df_pt <- df_pt %>%
  left_join(lengths_df_om, by = c("Assembly" = "Genome")) %>%
  rename(Reference_length = Length) %>%
  left_join(genome_map %>% select(Assembly, Species), by = "Assembly")

# Labels and filtering
label_map_pt <- c(
  "full_sim_pt_200000" = "NanoSim CM 200K",
  "full_sim_pt_100000" = "NanoSim CM 100K",
  "full_sim_pt_50000" = "NanoSim CM 50K",
  "full_sim_pt_20000" = "NanoSim CM 20K",
  "full_sim_pt0000" = "NanoSim CM 10K"
)

df_pt_long <- df_pt %>%
  pivot_longer(cols = all_of(assembler_cols_pt),
               names_to = "Dataset", values_to = "NGA50") %>%
  filter(Dataset %in% names(label_map_pt)) %>%
  mutate(Label = label_map_pt[Dataset]) %>%
  filter(Species %in% virus_species) %>%
  mutate(NGA50_Reference_Length = NGA50 / Reference_length)

df_pt_filtered <- df_pt_long %>%
  filter(Label %in% c("NanoSim CM 200K", "NanoSim CM 20K"))

# --- Combine OM and PT datasets ---
df_pt_filtered$Experiment <- "PT"
df_om_filtered$Experiment <- "OM"

df_combined <- bind_rows(df_pt_filtered,df_om_filtered)
df_combined <- df_combined %>%
  mutate(Experiment = recode(Experiment,
                             "PT" = "NanoSim Pretrained Model",
                             "OM" = "NanoSim Custom Model"))
df_combined <- df_combined %>%
  mutate(Assembler = factor(Assembler, levels = c("flye", "raven", "megahit")))
df_combined$Experiment <- factor(df_combined$Experiment,
                                 levels = c("NanoSim Pretrained Model", "NanoSim Custom Model"))


# --- Plot with facet strip --- "#1b9e77","#d95f02","#7570b3"
species_levels <- levels(factor(df_combined$Species))

y_positions <- seq_along(species_levels) - 0.5  # Midpoints between labels

ggplot(df_combined, aes(x = NGA50_Reference_Length, y = Species, color = Assembler, shape = Label)) +
  geom_point(position = position_dodge(width = 0.5), size = 4, alpha = 0.9, stroke = 0.4) +
  geom_hline(yintercept = y_positions, color = "grey80", linetype = "dashed", size = 0.3) +
  scale_x_continuous(breaks = seq(0, 1.0, by = 0.25), limits = c(0.000001, 1.0), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ Experiment, ncol = 1) +
  scale_color_manual(values = c("flye" = "#1B9E77", "raven" = "#D95F02", "megahit" = "#7570B3")) +
  scale_shape_manual(values = c("NanoSim CM 200K" = 16, "NanoSim CM 20K" = 17)) +
  theme_light(base_family = "Arial") +
  labs(
    x = "NGA50 / Reference Genome Length",
    y = NULL,
    color = "Assembler",
    shape = "Simulation Model"
  ) +
  theme(
    text = element_text(size = 11),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text = element_text(color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 10)
  )
