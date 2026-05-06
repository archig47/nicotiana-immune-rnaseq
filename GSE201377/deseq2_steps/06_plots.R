library(DESeq2)
library(ggplot2)
library(pheatmap)

dds <- readRDS("/users/fyp/fyp5/project/GSE201377/deseq2_steps/dds.rds")
res_shrunk <- readRDS("/users/fyp/fyp5/project/GSE201377/deseq2_steps/res_shrunk.rds")
outdir <- "/users/fyp/fyp5/project/GSE201377/deseq2"
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

vsd <- vst(dds, blind=TRUE)

# PCA
pca_data <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))
pdf(file.path(outdir, "PCA.pdf"), width=7, height=5)
ggplot(pca_data, aes(x=PC1, y=PC2, colour=condition, label=name)) +
  geom_point(size=4) + geom_text(vjust=-0.8, size=3) +
  xlab(paste0("PC1: ", pct_var[1], "%")) + ylab(paste0("PC2: ", pct_var[2], "%")) +
  ggtitle("PCA — GSE201377") + theme_bw()
dev.off()

# MA
pdf(file.path(outdir, "MA.pdf"), width=7, height=5)
plotMA(res_shrunk, main="MA plot — GSE201377 (apeglm)", ylim=c(-5,5))
dev.off()

# Volcano
res_df <- as.data.frame(res_shrunk)
res_df$sig <- ifelse(abs(res_df$log2FoldChange)>1 & res_df$padj<0.05,
                     ifelse(res_df$log2FoldChange>1, "Up in Psy", "Down in Psy"), "NS")
res_df_plot <- subset(res_df, !is.na(padj))
pdf(file.path(outdir, "volcano.pdf"), width=7, height=6)
ggplot(res_df_plot, aes(x=log2FoldChange, y=-log10(padj), colour=sig)) +
  geom_point(alpha=0.4, size=0.8) +
  scale_colour_manual(values=c("Up in Psy"="#e41a1c","Down in Psy"="#377eb8","NS"="grey60")) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  xlab("log2 Fold Change (Psy vs Pta)") + ylab("-log10(padj)") +
  ggtitle("Volcano — GSE201377") + theme_bw()
dev.off()

# Heatmap
sig <- subset(res_df, abs(log2FoldChange)>1 & padj<0.05)
top50 <- head(sig[order(sig$padj),], 50)
mat <- assay(vsd)[rownames(top50),]
mat <- mat - rowMeans(mat)
ann <- data.frame(condition=colData(dds)$condition, row.names=colnames(dds))
pdf(file.path(outdir, "heatmap_top50.pdf"), width=8, height=10)
pheatmap(mat, annotation_col=ann, fontsize_row=7, main="Top 50 DEGs — GSE201377")
dev.off()

cat("All plots saved to", outdir, "\n")
