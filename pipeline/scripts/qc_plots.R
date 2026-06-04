# QC summary plots for pipeline validation (Ma et al. 2025, PRJNA936199)
# Generates three bar plots from hardcoded MultiQC stats.
# To adapt for your own run: replace the values in the samples data frame below.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

out_dir <- "results/qc"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

condition_colours <- c("mock" = "#7F7F7F", "D36E" = "#C0392B", "dfliC" = "#2980B9")

samples <- data.frame(
  sample    = c("mock_rep1","mock_rep2","mock_rep3",
                "D36E_rep1","D36E_rep2","D36E_rep3",
                "dfliC_rep1","dfliC_rep2","dfliC_rep3"),
  condition = c("mock","mock","mock","D36E","D36E","D36E","dfliC","dfliC","dfliC"),
  read_M    = c(24.4, 24.7, 25.6, 23.4, 23.9, 26.2, 21.8, 20.2, 23.3),
  map_pct   = c(95.81, 94.95, 94.93, 95.44, 95.37, 95.51, 95.11, 94.75, 95.34),
  assigned  = c(41343373, 41506258, 42626166,
                39560250, 40927716, 44342371,
                36200873, 34057582, 38876601),
  stringsAsFactors = FALSE
)

# total mapped reads per sample for % assigned denominator
samples$total_reads <- samples$read_M * 1e6 * 2  # paired-end
samples$assigned_pct <- (samples$assigned / samples$total_reads) * 100

samples$condition <- factor(samples$condition, levels = c("mock","D36E","dfliC"))
samples$sample    <- factor(samples$sample, levels = samples$sample)

base_theme <- theme_classic(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    axis.title       = element_text(size = 10),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 11),
    legend.position  = "none",
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

# ── A: Read counts ─────────────────────────────────────────────────────────────
p_reads <- ggplot(samples, aes(x = sample, y = read_M, fill = condition)) +
  geom_col(width = 0.7, colour = "white") +
  geom_hline(yintercept = 20, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = condition_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Raw read depth", x = NULL, y = "Reads (millions)") +
  base_theme

# ── B: Unique mapping rate ─────────────────────────────────────────────────────
p_map <- ggplot(samples, aes(x = sample, y = map_pct, fill = condition)) +
  geom_col(width = 0.7, colour = "white") +
  geom_hline(yintercept = 75, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = condition_colours) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Uniquely mapped reads", x = NULL, y = "Uniquely mapped (%)") +
  base_theme

# ── C: % reads assigned to features ───────────────────────────────────────────
p_assigned <- ggplot(samples, aes(x = sample, y = assigned_pct, fill = condition)) +
  geom_col(width = 0.7, colour = "white") +
  geom_hline(yintercept = 70, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  scale_fill_manual(values = condition_colours) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Reads assigned to genes", x = NULL, y = "Assigned (%)") +
  base_theme

ggsave(file.path(out_dir, "qc_read_counts.pdf"),    p_reads,    width = 6, height = 4)
ggsave(file.path(out_dir, "qc_mapping_rate.pdf"),   p_map,      width = 6, height = 4)
ggsave(file.path(out_dir, "qc_assigned_reads.pdf"), p_assigned, width = 6, height = 4)

cat("QC plots saved to", out_dir, "\n")
