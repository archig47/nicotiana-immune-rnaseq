suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(ggrepel)
  library(yaml)
})

# ── Arguments ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
counts_file  <- args[1]   # results/counts/counts_matrix.tsv
samples_file <- args[2]   # samples.csv
config_file  <- args[3]   # config.yaml

config    <- yaml.load_file(config_file)
ref_cond  <- config$deseq2$reference_condition
lfc_thr   <- as.numeric(config$deseq2$lfc_threshold)
padj_thr  <- as.numeric(config$deseq2$padj_threshold)
top_n     <- as.integer(config$deseq2$top_n_heatmap)

out_dir   <- "results/deseq2"
plot_dir  <- file.path(out_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load counts ────────────────────────────────────────────────────────────────
raw <- read.table(counts_file, header = TRUE, skip = 1,
                  row.names = 1, check.names = FALSE)
# featureCounts columns: Chr Start End Strand Length [sample BAMs...]
counts <- raw[, 6:ncol(raw)]

# Clean column names: extract sample directory name (parent of the BAM file)
colnames(counts) <- basename(dirname(colnames(counts)))

# ── Load samples ───────────────────────────────────────────────────────────────
samples <- read.csv(samples_file, stringsAsFactors = FALSE)
samples <- samples[samples$sample_id %in% colnames(counts), ]

# Align count matrix columns to samples order
counts <- counts[, samples$sample_id]

col_data <- data.frame(
  row.names   = samples$sample_id,
  condition   = factor(samples$condition,
                       levels = c(ref_cond, setdiff(unique(samples$condition), ref_cond)))
)

# ── DESeq2 ─────────────────────────────────────────────────────────────────────
cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData   = col_data,
                               design    = ~ condition)

# Filter low-count genes: keep genes with >= 10 counts in at least min-group-size samples
min_samples <- min(table(col_data$condition))
keep <- rowSums(counts(dds) >= 10) >= min_samples
dds  <- dds[keep, ]
cat(sprintf("Genes after low-count filter: %d\n", sum(keep)))

dds <- DESeq(dds)

# ── Results (all contrasts vs reference) ──────────────────────────────────────
conditions <- levels(col_data$condition)
treatments <- conditions[conditions != ref_cond]

all_results <- list()
for (trt in treatments) {
  cat(sprintf("Extracting results: %s vs %s\n", trt, ref_cond))
  res <- results(dds, contrast = c("condition", trt, ref_cond))
  res <- lfcShrink(dds, contrast = c("condition", trt, ref_cond),
                   res = res, type = "ashr")
  res_df <- as.data.frame(res)
  res_df$gene_id   <- rownames(res_df)
  res_df$contrast  <- paste0(trt, "_vs_", ref_cond)
  res_df$direction <- ifelse(res_df$log2FoldChange > 0, "up", "down")
  all_results[[trt]] <- res_df
}

full_results <- do.call(rbind, all_results)
rownames(full_results) <- NULL

sig_results <- full_results[
  !is.na(full_results$padj) &
  full_results$padj < padj_thr &
  abs(full_results$log2FoldChange) > lfc_thr, ]

write.csv(full_results, file.path(out_dir, "DEGs_full.csv"),        row.names = FALSE)
write.csv(sig_results,  file.path(out_dir, "DEGs_significant.csv"), row.names = FALSE)
cat(sprintf("Significant DEGs: %d\n", nrow(sig_results)))

# ── Summary ────────────────────────────────────────────────────────────────────
sink(file.path(out_dir, "deseq2_summary.txt"))
cat("=== DESeq2 Summary ===\n\n")
cat(sprintf("Reference condition:   %s\n", ref_cond))
cat(sprintf("Contrasts run:         %s\n", paste(treatments, collapse=", ")))
cat(sprintf("Genes tested:          %d\n", nrow(dds)))
cat(sprintf("LFC threshold:         %.2f\n", lfc_thr))
cat(sprintf("padj threshold:        %.3f\n", padj_thr))
cat(sprintf("\nTotal significant DEGs: %d\n\n", nrow(sig_results)))
for (trt in treatments) {
  s <- sig_results[sig_results$contrast == paste0(trt, "_vs_", ref_cond), ]
  cat(sprintf("%s vs %s:\n", trt, ref_cond))
  cat(sprintf("  Upregulated:   %d\n", sum(s$log2FoldChange > 0)))
  cat(sprintf("  Downregulated: %d\n\n", sum(s$log2FoldChange < 0)))
}
sink()

# ── PCA plot ───────────────────────────────────────────────────────────────────
# vst() requires >= 1000 rows; fall back to varianceStabilizingTransformation for small datasets
vsd <- tryCatch(
  vst(dds, blind = TRUE),
  error = function(e) varianceStabilizingTransformation(dds, blind = TRUE)
)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                               colour = condition, label = name)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(size = 3, show.legend = FALSE) +
  labs(
    title    = "PCA of VST-normalised counts",
    x        = paste0("PC1: ", pct_var[1], "% variance"),
    y        = paste0("PC2: ", pct_var[2], "% variance"),
    colour   = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(plot_dir, "pca_plot.pdf"), p_pca, width = 7, height = 5)
cat("PCA plot saved.\n")

# ── Volcano plot (one per contrast) ───────────────────────────────────────────
for (trt in treatments) {
  res_df <- all_results[[trt]]
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df$neg_log10_p <- -log10(res_df$pvalue)
  res_df$neg_log10_p[!is.finite(res_df$neg_log10_p)] <- 0
  res_df$sig <- res_df$padj < padj_thr & abs(res_df$log2FoldChange) > lfc_thr
  res_df$category <- "NS"
  res_df$category[res_df$sig & res_df$log2FoldChange >  lfc_thr] <- "Up"
  res_df$category[res_df$sig & res_df$log2FoldChange < -lfc_thr] <- "Down"
  res_df$category <- factor(res_df$category, levels = c("NS", "Up", "Down"))

  top_labels <- res_df[res_df$sig, ]
  top_labels <- top_labels[order(top_labels$padj), ][seq_len(min(20, nrow(top_labels))), ]

  pal <- c("NS" = "grey80", "Up" = "#C0392B", "Down" = "#2980B9")
  n_up   <- sum(res_df$category == "Up")
  n_down <- sum(res_df$category == "Down")

  p_vol <- ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_p, colour = category)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = -log10(padj_thr), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr),    linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_text_repel(data = top_labels, aes(label = gene_id),
                    size = 2.5, max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = pal,
                        labels = c("NS", paste0("Up (", n_up, ")"), paste0("Down (", n_down, ")"))) +
    labs(
      title    = paste0("Volcano: ", trt, " vs ", ref_cond),
      subtitle = paste0("padj < ", padj_thr, "  |  |log2FC| > ", lfc_thr),
      x        = expression(log[2]~"Fold Change"),
      y        = expression(-log[10]~italic(p)),
      colour   = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(colour = "grey50", hjust = 0.5))

  fname <- if (length(treatments) == 1) "volcano_plot.pdf" else paste0("volcano_", trt, "_vs_", ref_cond, ".pdf")
  ggsave(file.path(plot_dir, fname), p_vol, width = 8, height = 6)
}
# Ensure the expected output path exists (for single-contrast case already done; multi-contrast: copy first)
if (length(treatments) > 1) {
  first_vol <- file.path(plot_dir, paste0("volcano_", treatments[1], "_vs_", ref_cond, ".pdf"))
  file.copy(first_vol, file.path(plot_dir, "volcano_plot.pdf"), overwrite = TRUE)
}
cat("Volcano plot(s) saved.\n")

# ── MA plot ────────────────────────────────────────────────────────────────────
for (trt in treatments) {
  res_df <- all_results[[trt]]
  res_df <- res_df[!is.na(res_df$padj), ]
  res_df$sig <- res_df$padj < padj_thr & abs(res_df$log2FoldChange) > lfc_thr

  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange,
                              colour = sig)) +
    geom_point(size = 0.6, alpha = 0.5) +
    geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    scale_colour_manual(values = c("FALSE" = "grey75", "TRUE" = "#E67E22"),
                        labels = c("Not significant", "DEG")) +
    labs(
      title  = paste0("MA plot: ", trt, " vs ", ref_cond),
      x      = expression(log[10]~"(mean expression + 1)"),
      y      = expression(log[2]~"Fold Change"),
      colour = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))

  fname <- if (length(treatments) == 1) "ma_plot.pdf" else paste0("ma_", trt, "_vs_", ref_cond, ".pdf")
  ggsave(file.path(plot_dir, fname), p_ma, width = 7, height = 5)
}
if (length(treatments) > 1) {
  first_ma <- file.path(plot_dir, paste0("ma_", treatments[1], "_vs_", ref_cond, ".pdf"))
  file.copy(first_ma, file.path(plot_dir, "ma_plot.pdf"), overwrite = TRUE)
}
cat("MA plot(s) saved.\n")

# ── Heatmap of top N DEGs ──────────────────────────────────────────────────────
top_genes <- sig_results[order(sig_results$padj), ]
top_genes <- top_genes[!duplicated(top_genes$gene_id), ]
top_genes <- head(top_genes, top_n)

if (nrow(top_genes) > 1) {
  mat <- assay(vsd)[top_genes$gene_id, , drop = FALSE]
  mat <- mat - rowMeans(mat)   # centre on row mean for visualisation

  ann_col <- data.frame(condition = col_data$condition,
                        row.names = rownames(col_data))

  pdf(file.path(plot_dir, "heatmap_top50.pdf"), width = 10, height = 12)
  pheatmap(mat,
           annotation_col  = ann_col,
           cluster_rows    = TRUE,
           cluster_cols    = TRUE,
           show_rownames   = nrow(mat) <= 80,
           fontsize_row    = 7,
           fontsize_col    = 9,
           color           = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
           main            = paste0("Top ", nrow(top_genes), " DEGs (VST, row-centred)"),
           border_color    = NA)
  dev.off()
  cat("Heatmap saved.\n")
} else {
  cat("Fewer than 2 significant DEGs — skipping heatmap.\n")
  pdf(file.path(plot_dir, "heatmap_top50.pdf"))
  plot.new()
  text(0.5, 0.5, "Not enough DEGs to plot heatmap", cex = 1.2)
  dev.off()
}

cat("\n=== DESeq2 analysis complete ===\n")
cat(paste0("Outputs written to: ", out_dir, "/\n"))
