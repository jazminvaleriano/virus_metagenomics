library(tidyverse)
library(ggalluvial)
library(stringr)

df <- read_tsv("04_positive_results/combined_report.tsv",col_select = c("Sample","Pass",
                                                                        "Lineage","Farm",
                                                                        "Symptoms","PCR"))%>%
                                                                        distinct()

# Add column for genus
df <- df %>%
  mutate(Genus = str_split(Lineage, ";", simplify = TRUE)[, 9])

df <- df %>%
  select(-Lineage) %>% distinct()

# List of viruses to exclude
excluded_taxa <- c(
  "Gemykibivirus anima1",
  "Gemykrogvirus carib1",
  "Porcine stool-associated circular virus/BEL/15V010",
  "Fur seal faeces associated circular DNA virus"
)

# Filter for valid detections, and exclude unwanted taxa
df_filtered <- df %>%
  filter(Pass == TRUE,
         !Genus %in% excluded_taxa)

df_filtered <- df_filtered %>%
  mutate(Genus = str_replace_all(Genus, c(
    "unclassified Mamastrovirus" = "Uncl. Mamastrovirus",
    "unclassified Parvovirinae" = "Uncl. Parvovirinae",
    "unclassified Porprismacovirus" = "Uncl. Porprismacovirus",
    "Porcine astrovirus 5" = "Porcine Astrovirus 5",
    "Orthocoronavirinae" = "Coronavirus")))

df_filtered$Genus <- fct_infreq(df_filtered$Genus)

# Group by farm, symptoms and Taxon
alluvial_data <- df_filtered %>%
  group_by(Farm, Symptoms, Genus) %>%
  summarise(Count = n(), .groups = "drop")


#library(ggthemes)  # Optional: for better themes

ggplot(alluvial_data,
       aes(axis1 = Farm, axis2 = Symptoms, axis3 = Genus, y = Count)) +
  scale_x_discrete(limits = c("Farm", "Symptoms", "Genus"),
                   expand = c(.05, .05)) +
  geom_alluvium(aes(fill = Genus), width = 1/5, alpha = 0.9) +
  geom_stratum(width = 1/8, color = "grey") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3.4, family = "sans", lineheight = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Viral Genus Detections - Farm → Symptoms → Genus",
    y = "Number of Detections",
    fill = "Genus"
  ) +
  guides(fill = guide_legend(ncol = 1))  # Stack legend vertically


#library(ggthemes)  # Optional: for better themes

ggplot(alluvial_data,
       aes(axis1 = Farm, axis2 = Symptoms, axis3 = Genus, y = Count)) +
  scale_x_discrete(limits = c("Farm", "Symptoms", "Genus"),
                   expand = c(.05, .05)) +
  geom_alluvium(aes(fill = Farm), width = 1/5, alpha = 0.9) +
  geom_stratum(width = 1/8, color = "grey") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3.4, family = "sans", lineheight = 0.8) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13, face = "bold"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Viral Genus Detections - Farm → Symptoms → Genus",
    y = "Number of Detections",
    fill = "Farm"
  ) +
  
  guides(fill = guide_legend(ncol = 1))  # Stack legend vertically

