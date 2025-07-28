# ─────────────────────────────────────────────
#         PERMANOVA FORWARD SELECTION
#        + SIDE-BY-SIDE COMPARISON PLOTS
# ─────────────────────────────────────────────

# Load libraries
library(vegan)
library(dplyr)
library(tibble)
library(ggplot2)
library(patchwork)

# ───── 1. LOAD & PREPARE DATA ─────

# Read metadata
meta <- read.csv("metadata.csv", stringsAsFactors = FALSE) %>%
  filter(Farm != "Control") %>%
  column_to_rownames("Sample") %>%
  select(Symptoms, Farm, Location, Age_bracket, Virus) %>%
  mutate(
    across(everything(), as.factor),
    Virus_Symptoms = interaction(Virus, Symptoms, drop = TRUE)
  )

# Load and clean beta diversity matrix
beta_mat <- read.table("02_diversity_analysis/beta_diversity_matrix.txt",
                       header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Remove control
beta_mat <- beta_mat[!rownames(beta_mat) %in% "barcode22", !colnames(beta_mat) %in% "barcode22"]
beta_dist <- as.dist(as.matrix(beta_mat))

# Match metadata to distance matrix
sample_names <- labels(beta_dist)
meta_subset <- meta[sample_names, ]

# ───── 2. DEFINE REUSABLE FUNCTION ─────
run_forward_selection <- function(beta_dist, meta_subset, current_vars = c(), 
                                  test_vars, plot_title = NULL, 
                                  suppress_x_axis = FALSE,
                                  suppress_y_axis = FALSE) 
  {
  results_list <- list()
  
  # Forward test loop
  for (var in test_vars) {
    trial_vars <- c(current_vars, var)
    trial_formula <- as.formula(paste("beta_dist ~", paste(trial_vars, collapse = " + ")))
    
    trial_model <- adonis2(trial_formula, data = meta_subset, permutations = 999, method = "bray")
    
    results_list[[length(results_list) + 1]] <- list(
      tested_variable = var,
      full_variable_set = trial_vars,
      R2 = trial_model$R2[1],
      p_value = trial_model$`Pr(>F)`[1],
      model = trial_model
    )
    
    cat("\n---- Testing:", paste(trial_vars, collapse = " + "), "----\n")
    print(trial_model)
  }
  
  # Create summary table
  summary_df <- do.call(rbind, lapply(results_list, function(x) {
    data.frame(
      tested_variable = x$tested_variable,
      model = paste(x$full_variable_set, collapse = " + "),
      R2 = x$R2,
      p_value = x$p_value
    )
  }))
  
  # Significance labels
  summary_df$Significance <- cut(summary_df$p_value,
                                 breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                                 labels = c("***", "**", "*", ".", "ns"))
  
  # Set x-axis label conditionally
  x_label <- if (suppress_x_axis) waiver() else "Explained Variation (R²)"
  
  # Create lollipop plot
  p <- ggplot(summary_df, aes(x = R2, y = reorder(model, R2), color = Significance)) +
    geom_segment(aes(x = 0, xend = R2, yend = model), size = 1.2) +
    geom_point(size = 5) +
    geom_text(aes(label = as.character(Significance)), hjust = -1.2, size = 5, na.rm = TRUE) +
    labs(
      title = plot_title,
      x = if (suppress_x_axis) NULL else "Explained Variation (R²)",
      y = if (suppress_y_axis) NULL else "Model"
    ) +
    xlim(0, 1) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.y = element_text(size = 12),
      plot.title = element_text(size = 16),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      legend.position = "none",
      axis.text.x = if (suppress_x_axis) element_blank() else element_text(),
      axis.ticks.x = if (suppress_x_axis) element_blank() else element_line(),
    ) +
    scale_color_manual(values = c("***" = "darkred", "**" = "red", "*" = "orange", "." = "gold", "ns" = "gray"),
                       na.value = "gray")
  
  return(p)
}


# ───── 3. RUN FORWARD SELECTION ROUNDS ─────

vars_all <- c("Farm", "Symptoms", "Location", "Age_bracket", "Virus", "Virus_Symptoms")
plot1 <- run_forward_selection(beta_dist, meta_subset,
                               current_vars = c(),
                               test_vars = vars_all,
                               plot_title = "Round 1: Starting from Null",
                               suppress_x_axis = TRUE,
                               suppress_y_axis = TRUE)

plot2 <- run_forward_selection(beta_dist, meta_subset,
                               current_vars = c("Farm"),
                               test_vars = setdiff(vars_all, "Farm"),
                               plot_title = "Round 2: Adding to Farm",
                               suppress_x_axis = TRUE,
                               suppress_y_axis = FALSE) # keep y axis here 

plot3 <- run_forward_selection(beta_dist, meta_subset,
                               current_vars = c("Farm", "Virus_Symptoms"),
                               test_vars = setdiff(vars_all, c("Farm", "Virus_Symptoms")),
                               plot_title = "Round 3: Adding to Farm + Virus_Symptoms",
                               suppress_x_axis = FALSE,
                               suppress_y_axis = TRUE)  # keep x-axis here


# ───── 4. COMBINE PLOTS SIDE BY SIDE ─────
combined_plot <- (plot1 / plot2 / plot3) +
  plot_layout(guides = "collect", ncol = 1) +
  plot_annotation(
    title = "PERMANOVA Forward Model Selection",
    caption = "Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1"
  ) &
  theme(
    plot.caption = element_text(size = 12, hjust = 0),
    legend.position = "none"
  )

combined_plot

# ───── 5. Save the figure ─────
ggsave("forward_selection_combined.png", combined_plot, width = 18, height = 6, dpi = 300)
