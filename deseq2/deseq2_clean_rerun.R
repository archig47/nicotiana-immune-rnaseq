suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

CLEAN_DIR  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
DESEQ_DIR  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2"
OUTLIER    <- "SRR23875062"

dir.create(file.path(CLEAN_DIR, "results"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CLEAN_DIR, "figures"),  showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CLEAN_DIR, "logs"),     showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(CLEAN_DIR, "nlr_prr"), showWarnings = FALSE, recursive = TRUE)

cat("=== Clean DESeq2 rerun — outlier removed ===\n")
cat(sprintf("Removing: %s\n\n", OUTLIER))

# ── 1. Load original counts and metadata ─────────────────────────────────
cat("Loading original dds...\n")
dds_old    <- readRDS(file.path(DESEQ_DIR, "dds_32.rds"))
counts_mat <- counts(dds_old)
cat(sprintf("Original: %d genes x %d samples\n", nrow(counts_mat), ncol(counts_mat)))

meta <- read_csv("/users/fyp/fyp5/project/PRJNA945175/clust/sample_metadata.csv",
                 show_col_types = FALSE)
meta <- meta[match(colnames(counts_mat), meta$sample), ]
stopifnot(all(meta$sample == colnames(counts_mat)))

# ── 2. Remove outlier ────────────────────────────────────────────────────
keep_samples  <- meta$sample != OUTLIER
meta_clean    <- meta[keep_samples, ]
counts_clean  <- counts_mat[, keep_samples]
cat(sprintf("After removal: %d genes x %d samples\n\n", nrow(counts_clean), ncol(counts_clean)))

# Save clean metadata
write_csv(meta_clean, file.path(CLEAN_DIR, "sample_metadata_clean.csv"))

# ── 3. Build clean DESeq2 object ─────────────────────────────────────────
cat("Building clean DESeq2 object...\n")
meta_clean$group <- factor(meta_clean$group)
meta_clean$group <- relevel(meta_clean$group, ref = "CK_1")

dds <- DESeqDataSetFromMatrix(
  countData = counts_clean,
  colData   = meta_clean,
  design    = ~ group
)

# Pre-filter
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat(sprintf("Genes after filtering: %d\n", nrow(dds)))

# ── 4. VST for PCA (before DESeq2 fit) ───────────────────────────────────
cat("Computing VST for PCA...\n")
vsd <- vst(dds, blind = TRUE)
saveRDS(vsd, file.path(CLEAN_DIR, "vsd_clean.rds"))

# Save VST matrix
vst_mat <- assay(vsd)
write.csv(vst_mat, file.path(CLEAN_DIR, "vst_matrix_clean.csv"))
cat("VST matrix saved\n")

# ── 5. PCA per timepoint ──────────────────────────────────────────────────
cat("Generating PCA plots...\n")

condition_colors <- c(
  CK   = "#888888",
  HCRV = "#2ecc71",
  TH   = "#9b59b6",
  TSWV = "#e74c3c"
)

make_pca_plot <- function(vsd_obj, timepoint_filter, label) {
  # Subset to timepoint
  tp_samples <- vsd_obj$timepoint == timepoint_filter
  vsd_sub    <- vsd_obj[, tp_samples]

  # PCA on top 500 variable genes
  rv      <- rowVars(assay(vsd_sub))
  top500  <- order(rv, decreasing = TRUE)[1:min(500, length(rv))]
  pca_res <- prcomp(t(assay(vsd_sub)[top500, ]))

  pct <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  pca_df <- data.frame(
    PC1       = pca_res$x[, 1],
    PC2       = pca_res$x[, 2],
    sample    = vsd_sub$sample,
    condition = vsd_sub$condition,
    group     = vsd_sub$group
  )

  ggplot(pca_df, aes(x = PC1, y = PC2, colour = condition)) +
    geom_point(size = 4, alpha = 0.9) +
    scale_colour_manual(values = condition_colors, name = "Condition") +
    labs(
      title    = label,
      x        = sprintf("PC1: %.1f%% variance", pct[1]),
      y        = sprintf("PC2: %.1f%% variance", pct[2])
    ) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(size = 11, colour = "grey40"),
      axis.text    = element_text(size = 9),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.3)
    )
}

p1 <- make_pca_plot(vsd, 1,  "1 dpi")
p2 <- make_pca_plot(vsd, 7,  "7 dpi")
p3 <- make_pca_plot(vsd, 14, "14 dpi")

# Combined figure with patchwork
library(patchwork)
p_combined <- (p1 + p2 + p3) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title    = "PCA of VST-transformed counts following outlier removal",
    subtitle = sprintf("Outlier removed: %s (CK 7 dpi, Cook's D = 6.30)", OUTLIER),
    theme = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey45")
    )
  )

ggsave(file.path(CLEAN_DIR, "figures", "pca_clean_per_timepoint.pdf"),
       p_combined, width = 12, height = 4.5, device = "pdf")
ggsave(file.path(CLEAN_DIR, "figures", "pca_clean_per_timepoint.png"),
       p_combined, width = 12, height = 4.5, dpi = 300)
cat("PCA saved\n")

# ── 6. Run DESeq2 ────────────────────────────────────────────────────────
cat("\nRunning DESeq2 (this may take a few minutes)...\n")
dds <- DESeq(dds)
saveRDS(dds, file.path(CLEAN_DIR, "dds_clean.rds"))
cat("dds_clean.rds saved\n")

# Save new size factors
sf_df <- data.frame(
  sample      = colnames(dds),
  size_factor = sizeFactors(dds)
)
write_csv(sf_df, file.path(CLEAN_DIR, "size_factors_clean.csv"))

# ── 7. Extract all 9 contrasts ───────────────────────────────────────────
conditions <- c("HCRV", "TSWV", "TH")
timepoints <- c(1, 7, 14)

extract_contrast <- function(dds, treatment, tp) {
  treat_group <- paste0(treatment, "_", tp)
  ck_group    <- paste0("CK_", tp)
  cat(sprintf("  %s vs %s\n", treat_group, ck_group))

  res <- results(dds,
                 contrast = c("group", treat_group, ck_group),
                 alpha    = 0.05)

  as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    arrange(padj)
}

cat("\nExtracting timepoint-specific contrasts...\n")
deg_summary <- data.frame()

for (cond in conditions) {
  for (tp in timepoints) {
    res_df <- extract_contrast(dds, cond, tp)

    out_file <- file.path(CLEAN_DIR, "results",
                          sprintf("results_%s_%sdpi_vs_CK.csv", cond, tp))
    write_csv(res_df, out_file)

    n_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange >  1, na.rm = TRUE)
    n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
    n_deg  <- n_up + n_down

    cat(sprintf("    DEGs: %d (%d up, %d down)\n", n_deg, n_up, n_down))

    deg_summary <- rbind(deg_summary, data.frame(
      contrast  = sprintf("%s_%sdpi_vs_CK", cond, tp),
      virus     = cond,
      timepoint = tp,
      n_up      = n_up,
      n_down    = n_down,
      n_deg     = n_deg
    ))
  }
}

write_csv(deg_summary, file.path(CLEAN_DIR, "deg_summary.csv"))
cat("\nDEG summary:\n")
print(deg_summary)

cat("\n=== Done. All outputs in deseq2_clean/ ===\n")
