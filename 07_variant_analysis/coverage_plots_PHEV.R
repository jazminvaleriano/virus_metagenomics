library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rtracklayer)


## For many bams
bam_files <- c("02_medaka/barcode33_Porcine_Hemagglutinating/calls_to_ref.bam", 
               "02_medaka/barcode34_Porcine_Hemagglutinating/calls_to_ref.bam", 
               "02_medaka/barcode35_Porcine_Hemagglutinating/calls_to_ref.bam")
sample_names <- c("Barcode 33", "Barcode 34", "Barcode 35")
gff_file <- "Hemagglutinating_CoV/Porcine_Hemagglutinating_genome/sequence.gff3"


contig <- "OP959790.1"  # your reference name

# Loop over BAMs and get coverage
all_cov <- lapply(seq_along(bam_files), function(i) {
  bam <- BamFile(bam_files[i])
  cov <- coverage(bam)[[contig]]
  data.frame(
    position = seq_along(cov),
    coverage = as.numeric(cov),
    sample = sample_names[i]
  )
})

# Combine into a single dataframe
cov_df <- bind_rows(all_cov)

# ───── get features from Gff ─────
gff <- import(gff_file)
cds <- gff[gff$type == "CDS"]
cds_df <- as.data.frame(cds)

# ───── Group by gene (some are repeated) ─────
annotations <- cds_df %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise(
    start = min(start),
    end = max(end),
    strand = unique(strand),
    .groups = "drop"
  )
# ───── Avoid genes overlap ─────
annotations <- annotations[order(annotations$start), ]
annotations$track <- NA_integer_

for (i in seq_len(nrow(annotations))) {
  current <- annotations[i, ]
  for (track in 1:100) {
    overlapping <- which(annotations$track == track &
                           !(annotations$end < current$start | annotations$start > current$end))
    if (length(overlapping) == 0) {
      annotations$track[i] <- track
      break
    }
  }
}

# ───── Set vertical positions ─────
annotations$ymin <- -5 - (annotations$track * 5)
annotations$ymax <- annotations$ymin + 3

# ───── Gene labels (manual positioning for short ones) ─────
annotations$label <- annotations$gene
annotations$label_y <- (annotations$ymin + annotations$ymax) / 2  # default center

# Move some genes label above or below
annotations$label_y[annotations$gene == "NS4.9"] <- annotations$ymax[annotations$gene == "NS4.9"] + 1
annotations$label_y[annotations$gene == "NS12.7"] <- annotations$ymax[annotations$gene == "NS12.7"] + 1
annotations$label_y[annotations$gene == "NS2"] <- annotations$ymax[annotations$gene == "NS2"] + 1

# ───── Coverage plot ─────
# Add log-transformed coverage (with pseudocount to handle zeros)
cov_df$log_coverage <- log10(cov_df$coverage + 1)

# Define coverage plot (log scale)
coverage_log <- ggplot(cov_df, aes(x = position, y = log_coverage)) +
  geom_area(fill = "#D6B7C9", alpha = 0.8) +
  geom_line(color = "#3B1628", linewidth = 0.3) +    
  facet_wrap(~sample, ncol = 1, scales = "free_x") +
  scale_y_continuous(
    name = "Coverage",
    breaks = log10(c(1, 2, 5, 10, 20, 50, 100) + 1),
    labels = c("1", "2", "5", "10", "20", "50", "100"),
    expand = expansion(mult = c(0.01, 0.1)),  # add some breathing room
    minor_breaks = NULL
  ) +
  theme_minimal(base_size = 14) +
  labs(x = NULL) +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
  )

# ───── Annotation plot ─────
annotation_plot <- ggplot() +
  geom_rect(data = annotations,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene),
            color = "black", linewidth = 0.2, alpha = 1) +
  geom_text(data = annotations,
            aes(x = (start + end)/2, y = label_y, label = label),
            size = 5, # this is the text size
            hjust = 0.5, angle = 0, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)  #  Add space to top of genes for labels
    
  ) +
  coord_cartesian(ylim = c(min(annotations$ymin), max(annotations$ymax) + 1.5))

# ───── Combine plots ─────
final_plot <- (coverage_log / annotation_plot) +
  plot_layout(heights = c(3, 0.7))

final_plot

