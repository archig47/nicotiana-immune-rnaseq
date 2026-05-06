library(DESeq2)
library(apeglm)

dds <- readRDS("/users/fyp/fyp5/project/GSE201377/deseq2_steps/dds.rds")

res_shrunk <- lfcShrink(dds, coef="condition_Psy_vs_Pta", type="apeglm")
summary(res_shrunk)

res_df <- as.data.frame(res_shrunk)
sig <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05)
cat("Total DEGs:", nrow(sig), "\n")
cat("Up in Psy:", nrow(subset(sig, log2FoldChange > 1)), "\n")
cat("Down in Psy:", nrow(subset(sig, log2FoldChange < -1)), "\n")

saveRDS(res_shrunk, "/users/fyp/fyp5/project/GSE201377/deseq2_steps/res_shrunk.rds")
cat("res_shrunk saved.\n")
