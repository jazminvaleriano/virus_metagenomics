library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(ggplot2)
library(patchwork)
library(rtracklayer)

# ───── Input files ─────
bam_files <- c(
  "02_medaka/barcode35_swIAV_H1N1_FELDBACH/calls_to_ref.bam",
  "02_medaka/barcode36_swIAV_H1N1_FELDBACH/calls_to_draft.bam",
  "02_medaka/barcode47_swIAV_H1N1_FELDBACH/calls_to_draft.bam"
  
)
sample_names <- c("Barcode 35", "Barcode 36","Barcode 47")
gff_file <- "genomes/H1N1_Feldbach/sequence (1).gff3"

# ───── GFF import and preprocessing ─────
gff <- import(gff_file) %>% as.data.frame()
gff$contig <- as.character(gff$seqnames)

# ───── BAM coverage for all contigs ─────
all_cov <- lapply(seq_along(bam_files), function(i) {
  bam <- BamFile(bam_files[i])
  cov_list <- coverage(bam)
  do.call(rbind, lapply(names(cov_list), function(ctg) {
    cov <- cov_list[[ctg]]
    data.frame(
      position = seq_along(cov),
      coverage = as.numeric(cov),
      sample = sample_names[i],
      contig = ctg
    )
  }))
})
cov_df <- bind_rows(all_cov)

# ───── Calculate contig offsets for global positioning ─────
contig_lengths <- cov_df %>%
  group_by(contig) %>%
  summarise(length = max(position), .groups = "drop") %>%
  mutate(x_offset = lag(cumsum(length), default = 0))

# ───── Add vertical dashed lines at contig boundaries ─────
boundary_lines <- contig_lengths %>%
  filter(x_offset > 0)  # don't draw a line at 0


# ───── Apply global x-coordinates ─────
cov_df <- cov_df %>%
  left_join(contig_lengths, by = "contig") %>%
  mutate(global_pos = position + x_offset,
         log_coverage = log10(coverage + 1))

# ───── Annotations: CDS grouped by gene ─────
annotations <- gff %>%
  filter(type == "CDS", !is.na(gene)) %>%
  group_by(contig, gene) %>%
  summarise(
    start = min(start),
    end = max(end),
    strand = unique(strand),
    .groups = "drop"
  ) %>%
  left_join(contig_lengths, by = "contig") %>%
  mutate(
    start = start + x_offset,
    end = end + x_offset
  )

# ───── Assign non-overlapping tracks ─────
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

# ───── Vertical positions and labels ─────
annotations$ymin <- -5 - (annotations$track * 5)
annotations$ymax <- annotations$ymin + 3
annotations$label <- annotations$gene
annotations$label_y <- (annotations$ymin + annotations$ymax) / 2

# Optional manual label shifts
annotations$label_y[annotations$gene == "NS4.9"] <- annotations$ymax[annotations$gene == "NS4.9"] + 1
annotations$label_y[annotations$gene == "NS12.7"] <- annotations$ymax[annotations$gene == "NS12.7"] + 1
annotations$label_y[annotations$gene == "NS2"] <- annotations$ymax[annotations$gene == "NS2"] + 1
# ───── Contig (accession) axis labels ─────
contig_labels <- contig_lengths %>%
  mutate(center = x_offset + length / 2)


# ───── Plot: Coverage ─────
coverage_plot <- ggplot(cov_df, aes(x = global_pos, y = log_coverage, fill = sample)) +
  geom_area(alpha = 0.3, position = "identity") +
  geom_line(aes(color = sample), linewidth = 0.3) +
  geom_vline(data = boundary_lines, aes(xintercept = x_offset),
             linetype = "dashed", color = "grey30", linewidth = 0.4) +
  scale_y_continuous(
    name = "Coverage",
    breaks = log10(c(1, 2, 5, 10, 20, 50, 100) + 1),
    labels = c("1", "2", "5", "10", "20", "50", "100"),
    expand = expansion(mult = c(0.01, 0.1)),
    minor_breaks = NULL
  ) +
  scale_x_continuous(
    breaks = contig_labels$center,
    labels = contig_labels$contig
  ) +
  theme_minimal(base_size = 14) +
  labs(x = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# ───── Plot: Annotations ─────
annotation_plot <- ggplot() +
  geom_rect(data = annotations,
            aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = gene),
            color = "black", linewidth = 0.2, alpha = 1) +
  geom_text(data = annotations,
            aes(x = (start + end)/2, y = label_y, label = label),
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
  coord_cartesian(ylim = c(min(annotations$ymin), max(annotations$ymax) + 1.5))

coverage_plot <- coverage_plot +
  scale_x_continuous(
    breaks = contig_labels$center,
    labels = contig_labels$contig
  )

# ───── Combine plots ─────
final_plot <- (coverage_plot / annotation_plot) +
  plot_layout(heights = c(3, 0.7))

# ───── Show plot ─────
final_plot
