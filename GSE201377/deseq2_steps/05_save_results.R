library(DESeq2)

res_shrunk <- readRDS("/users/fyp/fyp5/project/GSE201377/deseq2_steps/res_shrunk.rds")
outdir <- "/users/fyp/fyp5/project/GSE201377/deseq2"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

res_df <- as.data.frame(res_shrunk)
res_df$gene_id <- rownames(res_df)
sig <- subset(res_df, abs(log2FoldChange) > 1 & padj < 0.05)

write.csv(res_df, file.path(outdir, "all_results.csv"),        row.names=FALSE)
write.csv(sig,    file.path(outdir, "significant_DEGs.csv"),   row.names=FALSE)
write.csv(subset(sig, log2FoldChange >  1), file.path(outdir, "upregulated_Psy.csv"),   row.names=FALSE)
write.csv(subset(sig, log2FoldChange < -1), file.path(outdir, "downregulated_Psy.csv"), row.names=FALSE)

cat("Files written to", outdir, "\n")
