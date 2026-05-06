library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

setwd("/users/fyp/fyp5/project/PRJNA945175/deseq2")
dir.create("plots", showWarnings = FALSE)

dds <- readRDS("dds_32.rds")
rld <- readRDS("rld_32.rds")

# Clean sample names
clean_names <- function(x) {
  x <- gsub(".*SRR", "SRR", x)
  x <- gsub("\\..*", "", x)
  x
}
colnames(rld) <- clean_names(colnames(rld))
colnames(dds) <- clean_names(colnames(dds))

# Load results safely
read_res <- function(f) {
  df <- read.csv(f)
  rownames(df) <- make.unique(as.character(df$gene_id))
  df
}
res_HCRV <- read_res("all_results32_HCRV_vs_CK.csv")
res_TSWV <- read_res("all_results32_TSWV_vs_CK.csv")
res_TH   <- read_res("all_results32_TH_vs_CK.csv")

nlr_HCRV <- read_res("DEGs32_NLR_HCRV_vs_CK.csv")
nlr_TSWV <- read_res("DEGs32_NLR_TSWV_vs_CK.csv")
nlr_TH   <- read_res("DEGs32_NLR_TH_vs_CK.csv")

prr_HCRV <- read_res("DEGs32_PRR_HCRV_vs_CK.csv")
prr_TSWV <- read_res("DEGs32_PRR_TSWV_vs_CK.csv")
prr_TH   <- read_res("DEGs32_PRR_TH_vs_CK.csv")

# 1. PCA
pca_data <- plotPCA(rld, intgroup=c("condition","timepoint"), returnData=TRUE)
pct <- round(100 * attr(pca_data, "percentVar"))

pdf("plots/PCA.pdf", width=8, height=6)
print(ggplot(pca_data, aes(PC1, PC2, color=condition, shape=timepoint)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", pct[1], "% variance")) +
  ylab(paste0("PC2: ", pct[2], "% variance")) +
  ggtitle("PCA - PRJNA945175") +
  theme_bw())
dev.off()

# 2. Sample correlation heatmap
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
ann <- as.data.frame(colData(dds)[, c("condition","timepoint")])

pdf("plots/sample_correlation_heatmap.pdf", width=10, height=9)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_col=ann,
         col=colors,
         main="Sample distances - PRJNA945175")
dev.off()

# 3. Library sizes
lib_df <- data.frame(
  sample=colnames(dds),
  mapped=colSums(counts(dds)),
  condition=colData(dds)$condition
)

pdf("plots/library_sizes.pdf", width=10, height=5)
print(ggplot(lib_df, aes(x=reorder(sample, mapped), y=mapped/1e6, fill=condition)) +
  geom_bar(stat="identity") +
  coord_flip() +
  xlab("Sample") +
  ylab("Mapped reads (millions)") +
  ggtitle("Library sizes - PRJNA945175") +
  theme_bw())
dev.off()

# 4. Dispersion plot
pdf("plots/dispersion_plot.pdf", width=7, height=5)
plotDispEsts(dds, main="Dispersion estimates - PRJNA945175")
dev.off()

# 5. Volcano plots
make_volcano <- function(res, title, filename) {
  df <- as.data.frame(res)
  df$gene <- rownames(df)
  df$DEG <- ifelse(!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1, "DEG", "NS")
  top <- df[df$DEG == "DEG",]
  if (nrow(top) > 0) {
    top <- head(top[order(top$padj),], 10)
  }
  pdf(filename, width=8, height=6)
  p <- ggplot(df, aes(log2FoldChange, -log10(padj), color=DEG)) +
    geom_point(alpha=0.5, size=0.8) +
    scale_color_manual(values=c(DEG="red", NS="grey60")) +
    geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey40") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
    ggtitle(title) +
    theme_bw()
  if (nrow(top) > 0) {
    p <- p + geom_text(data=top, aes(label=gene), size=2.5, vjust=-0.5, color="black")
  }
  print(p)
  dev.off()
}

make_volcano(res_HCRV, "Volcano - HCRV vs CK", "plots/volcano_HCRV.pdf")
make_volcano(res_TSWV, "Volcano - TSWV vs CK", "plots/volcano_TSWV.pdf")
make_volcano(res_TH,   "Volcano - TH vs CK",   "plots/volcano_TH.pdf")

# 6. Heatmap top 50 DEGs
get_top_degs <- function(res, n=20) {
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj) & df$padj < 0.05,]
  head(rownames(df[order(df$padj),]), n)
}

degs_all <- unique(c(
  get_top_degs(res_HCRV),
  get_top_degs(res_TSWV),
  get_top_degs(res_TH)
))
degs_all <- degs_all[degs_all %in% rownames(assay(rld))]
degs_all <- degs_all[1:min(50, length(degs_all))]

mat <- assay(rld)[degs_all,]
mat <- t(scale(t(mat)))
ann <- as.data.frame(colData(dds)[, c("condition","timepoint")])

pdf("plots/heatmap_top50_DEGs.pdf", width=12, height=10)
pheatmap(mat,
         annotation_col=ann,
         show_rownames=TRUE,
         show_colnames=FALSE,
         main="Top 50 DEGs - PRJNA945175",
         fontsize_row=6)
dev.off()

# 7. NLR dotplots
make_dotplot <- function(res_sub, title, filename, color_high="firebrick") {
  df <- as.data.frame(res_sub)
  df$gene <- rownames(df)
  df <- df[order(df$baseMean, decreasing=TRUE),]
  df <- head(df, 20)
  df$absLFC <- abs(df$log2FoldChange)
  pdf(filename, width=7, height=7)
  p <- ggplot(df, aes(x=log2FoldChange, y=reorder(gene, log2FoldChange),
                      size=log10(baseMean), color=absLFC)) +
    geom_point() +
    scale_color_gradient(low="grey80", high=color_high) +
    geom_vline(xintercept=0, linetype="dashed") +
    xlab("log2 Fold Change") +
    ylab("Gene ID") +
    ggtitle(title) +
    theme_bw()
  print(p)
  dev.off()
}

make_dotplot(nlr_HCRV, "Top NLRs - HCRV vs CK", "plots/NLR_dotplot_HCRV.pdf")
make_dotplot(nlr_TSWV, "Top NLRs - TSWV vs CK", "plots/NLR_dotplot_TSWV.pdf")
make_dotplot(nlr_TH,   "Top NLRs - TH vs CK",   "plots/NLR_dotplot_TH.pdf")

make_dotplot(prr_HCRV, "Top PRRs - HCRV vs CK", "plots/PRR_dotplot_HCRV.pdf", color_high="steelblue")
make_dotplot(prr_TSWV, "Top PRRs - TSWV vs CK", "plots/PRR_dotplot_TSWV.pdf", color_high="steelblue")
make_dotplot(prr_TH,   "Top PRRs - TH vs CK",   "plots/PRR_dotplot_TH.pdf",   color_high="steelblue")

message("All plots saved to plots/")
