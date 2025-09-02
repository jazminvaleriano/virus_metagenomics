library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rtracklayer)

# ───── Input files ─────
bam_files <- c(
  "02_medaka/barcode25_Porcine_Hokovirus//calls_to_ref.bam",
  "02_medaka/barcode26_Porcine_Hokovirus//calls_to_ref.bam"
)
sample_names <- c("Barcode 25", "Barcode 26")
gff_file <- "Parvoviridae/genomes/Porcine_Hokovirus/Porcine_Hokovirus.gff3"
contig <- "NC_038546.1"  # Your single accession name

# ───── Load GFF and extract CDS annotations ─────
gff <- import(gff_file)
gff_df <- as.data.frame(gff)
gff_df <- gff_df %>%
  filter(seqnames == contig, type == "CDS", !is.na(gene)) %>%
  group_by(gene) %>%
  summarise(
    start = min(start, na.rm = TRUE),
    end = max(end, na.rm = TRUE),
    strand = unique(strand),
    .groups = "drop"
  )

# ───── Assign non-overlapping tracks (for stackinggene# ───── Assign non-overlapping tracks (for stacking) ─────
gff_df <- gff_df[order(gff_df$start), ]
gff_df$track <- NA_integer_

for (i in seq_len(nrow(gff_df))) {
  current <- gff_df[i, ]
  for (track in 1:100) {
    overlapping <- which(gff_df$track == track &
                           !(gff_df$end < current$start | gff_df$start > current$end))
    if (length(overlapping) == 0) {
      gff_df$track[i] <- track
      break
    }
  }
}

gff_df$ymin <- -5 - (gff_df$track * 5)
gff_df$ymax <- gff_df$ymin + 3
gff_df$label_y <- (gff_df$ymin + gff_df$ymax) / 2

# ───── Load BAM coverage ─────
all_cov <- lapply(seq_along(bam_files), function(i) {
  bam <- BamFile(bam_files[i])
  cov <- coverage(bam)[[contig]]
  data.frame(
    position = seq_along(cov),
    coverage = as.numeric(cov),
    sample = sample_names[i]
  )
})
cov_df <- bind_rows(all_cov)
cov_df$log_coverage <- log10(cov_df$coverage + 1)

# ───── Plot coverage ─────
coverage_plot <- ggplot(cov_df, aes(x = position, y = log_coverage, fill = sample)) +
  geom_area(alpha = 0.4, position = "identity") +
  geom_line(aes(color = sample), linewidth = 0.3) +
  scale_y_continuous(
    name = "Coverage",
    breaks = log10(c(1, 2, 5, 10, 20, 50, 100) + 1),
    labels = c("1", "2", "5", "10", "20", "50", "100"),
    expand = expansion(mult = c(0.01, 0.1))
  ) +
  theme_minimal(base_size = 14) +
  labs(x = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ───── Plot annotations ─────
annotation_plot <- ggplot() +
  geom_rect(data = gff_df,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene),
            color = "black", linewidth = 0.2) +
  geom_text(data = gff_df,
            aes(x = (start + end)/2, y = label_y, label = gene),
            size = 5, hjust = 0.5, color = "black") +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal(base_size = 14) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),
    plot.margin = margin(t = 20, r = 5, b = 5, l = 5)
  ) +
  coord_cartesian(ylim = c(min(gff_df$ymin), max(gff_df$ymax) + 1.5))

# ───── Combine plots ─────
final_plot <- (coverage_plot / annotation_plot) +
  plot_layout(heights = c(3, 0.7))

# ───── Show plot ─────
final_plot
