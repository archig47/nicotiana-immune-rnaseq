library(DESeq2)
library(apeglm)

counts_raw <- read.table(
  "/users/fyp/fyp5/project/GSE201377/featurecounts_softmasked/counts_matrix_v12_multimapper.txt",
  header = TRUE, skip = 1, sep = "\t", row.names = 1
)

counts <- counts_raw[, 6:ncol(counts_raw)]
colnames(counts) <- gsub(".*/GSM[0-9]+_([^/]+)/Aligned.*", "\\1", colnames(counts))
rownames(counts) <- gsub("-mRNA$", "", rownames(counts))
counts <- round(counts)

cat("Dimensions before filtering:", dim(counts), "\n")
keep <- rowSums(counts) >= 10
counts <- counts[keep, ]
cat("Dimensions after filtering (>=10 counts):", dim(counts), "\n")

coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c(
    "Pta","Pta","Pta",
    "Psy","Psy","Psy",
    "Pta","Pta","Pta",
    "Psy","Psy","Psy"
  ), levels = c("Pta", "Psy"))
)

cat("Sample metadata:\n")
print(coldata)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

res_shrunk <- lfcShrink(dds, coef = "condition_Psy_vs_Pta", type = "apeglm")

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)
res_ordered <- res_df[order(res_df$padj), ]

degs      <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05)
degs_up   <- subset(degs, log2FoldChange > 1)
degs_down <- subset(degs, log2FoldChange < -1)

cat("\n── DEG Summary ──────────────────────────────────────────────────────\n")
cat("Total DEGs (|LFC|>1, padj<0.05):", nrow(degs), "\n")
cat("Upregulated in Psy:             ", nrow(degs_up), "\n")
cat("Downregulated in Psy:           ", nrow(degs_down), "\n")

outdir <- "/users/fyp/fyp5/project/GSE201377/deseq2_softmasked"
write.csv(res_ordered, file.path(outdir, "all_results.csv"), row.names = FALSE)
write.csv(degs,        file.path(outdir, "DEGs_sig.csv"),    row.names = FALSE)
write.csv(degs_up,     file.path(outdir, "DEGs_up_Psy.csv"), row.names = FALSE)
write.csv(degs_down,   file.path(outdir, "DEGs_down_Psy.csv"), row.names = FALSE)
saveRDS(dds,        file.path(outdir, "dds.rds"))
saveRDS(res_shrunk, file.path(outdir, "res_shrunk.rds"))

cat("\nAll outputs saved to", outdir, "\n")
