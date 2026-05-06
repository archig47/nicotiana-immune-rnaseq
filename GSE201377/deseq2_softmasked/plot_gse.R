library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

outdir <- "/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/plots"
dir.create(outdir, showWarnings=FALSE)

# ── Load objects ──────────────────────────────────────────────────────────
cat("Loading DESeq2 objects...\n")
dds       <- readRDS("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/dds.rds")
res       <- readRDS("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/res_shrunk.rds")
degs      <- read.csv("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_annotated.csv")
nlr_degs  <- read.csv("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_NLR.csv")
all_res   <- read.csv("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/all_results.csv")

# NLRtracker for expression dotplot
nlr <- read.table("/users/fyp/fyp5/project/genome/jiorgos/NbT2T.NLRtracker_classified.tsv",
                  header=TRUE, sep="\t")
nlr$gene_id <- gsub("-mRNA$", "", nlr$protein_id)

# rlog transform
rld <- rlog(dds, blind=FALSE)

# Condition colours
cond_cols <- c(Pta="#2166ac", Psy="#d6604d")
condition <- colData(dds)$condition

# ── 1. PCA plot ───────────────────────────────────────────────────────────
cat("PCA plot...\n")
pca_data <- plotPCA(rld, intgroup="condition", returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color=condition, label=name)) +
  geom_point(size=4) +
  ggrepel::geom_text_repel(size=3, max.overlaps=20) +
  scale_color_manual(values=cond_cols) +
  labs(title="PCA — GSE201377 (Psy vs Pta)",
       x=paste0("PC1: ", pct_var[1], "% variance"),
       y=paste0("PC2: ", pct_var[2], "% variance")) +
  theme_bw(base_size=13)
ggsave(file.path(outdir, "PCA.pdf"), p, width=7, height=5)

# ── 2. Sample-to-sample correlation heatmap ───────────────────────────────
cat("Correlation heatmap...\n")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
anno_col <- data.frame(condition=colData(rld)$condition,
                       row.names=colnames(rld))
anno_colors <- list(condition=cond_cols)

pdf(file.path(outdir, "sample_correlation_heatmap.pdf"), width=8, height=7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         annotation_col=anno_col,
         annotation_colors=anno_colors,
         main="Sample-to-sample distances — GSE201377")
dev.off()

# ── 3. Dispersion plot ────────────────────────────────────────────────────
cat("Dispersion plot...\n")
pdf(file.path(outdir, "dispersion_plot.pdf"), width=7, height=5)
plotDispEsts(dds, main="Dispersion estimates — GSE201377")
dev.off()

# ── 4. Volcano plot ───────────────────────────────────────────────────────
cat("Volcano plot...\n")
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                     "DEG", "NS")

# Key genes to label
key_genes <- c("NbCDPK2-1","NbCDPK2-2","NbRIPK-2","NbSRG3-1","NbSRG3-2",
                "NbLDOX-2","NbERF6-3","NbMBF1C")
res_df$label <- ifelse(res_df$gene_id %in% degs$gene_id[degs$gene_name %in% key_genes],
                        degs$gene_name[match(res_df$gene_id, degs$gene_id)], NA)
# Also label ATP synthase
res_df$label[res_df$gene_id == "NbT2T12g01957"] <- "ATP synthase α"

p <- ggplot(res_df, aes(log2FoldChange, -log10(padj), color=sig)) +
  geom_point(alpha=0.5, size=0.8) +
  scale_color_manual(values=c(DEG="#d6604d", NS="grey70")) +
  ggrepel::geom_text_repel(aes(label=label), size=3, max.overlaps=20,
                            na.rm=TRUE, color="black") +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  labs(title="Volcano plot — Psy vs Pta (GSE201377)",
       x="log2 Fold Change (Psy/Pta)",
       y="-log10(adjusted p-value)",
       color="") +
  theme_bw(base_size=13)
ggsave(file.path(outdir, "volcano.pdf"), p, width=8, height=6)

# ── 5. MA plot ────────────────────────────────────────────────────────────
cat("MA plot...\n")
pdf(file.path(outdir, "MA_plot.pdf"), width=7, height=5)
plotMA(res, main="MA plot — Psy vs Pta (GSE201377)", ylim=c(-8,8))
dev.off()

# ── 6. Heatmap top 50 DEGs ────────────────────────────────────────────────
cat("Top DEG heatmap...\n")
top50 <- degs %>% arrange(padj) %>% head(50) %>% pull(gene_id)
mat   <- assay(rld)[top50, ]
mat   <- mat - rowMeans(mat)  # Z-score
rownames(mat) <- degs$gene_name[match(rownames(mat), degs$gene_id)]
rownames(mat)[is.na(rownames(mat))] <- top50[is.na(rownames(mat))]

anno_col <- data.frame(condition=colData(rld)$condition, row.names=colnames(rld))
pdf(file.path(outdir, "heatmap_top50_DEGs.pdf"), width=9, height=12)
pheatmap(mat,
         annotation_col=anno_col,
         annotation_colors=anno_colors,
         show_colnames=TRUE,
         fontsize_row=8,
         scale="none",
         main="Top 50 DEGs — Psy vs Pta (GSE201377)")
dev.off()

# ── 7. p-value histogram ──────────────────────────────────────────────────
cat("p-value histogram...\n")
p <- ggplot(res_df %>% filter(!is.na(pvalue)), aes(pvalue)) +
  geom_histogram(bins=50, fill="#2166ac", color="white") +
  labs(title="p-value distribution — GSE201377",
       x="p-value", y="Count") +
  theme_bw(base_size=13)
ggsave(file.path(outdir, "pvalue_histogram.pdf"), p, width=6, height=4)

# ── 8. NLR expression dotplot ─────────────────────────────────────────────
cat("NLR expression dotplot...\n")
# Get top 20 NLRs by baseMean from all_results
nlr_expr <- all_res %>%
  filter(gene_id %in% nlr$gene_id) %>%
  arrange(desc(baseMean)) %>%
  head(20)

p <- ggplot(nlr_expr, aes(x=log2FoldChange, y=reorder(gene_id, log2FoldChange),
                           color=log2FoldChange, size=log10(baseMean+1))) +
  geom_point() +
  scale_color_gradient2(low="#2166ac", mid="grey80", high="#d6604d", midpoint=0) +
  geom_vline(xintercept=0, linetype="dashed", color="grey40") +
  labs(title="Top 20 NLRs by expression — GSE201377",
       x="log2 Fold Change (Psy/Pta)",
       y="NLR gene ID",
       size="log10(baseMean)",
       color="LFC") +
  theme_bw(base_size=11)
ggsave(file.path(outdir, "NLR_expression_dotplot.pdf"), p, width=8, height=7)

# ── 9. Library size bar plot ──────────────────────────────────────────────
cat("Library size bar plot...\n")
lib_sizes <- colSums(counts(dds))
lib_df <- data.frame(
  sample    = names(lib_sizes),
  counts    = lib_sizes,
  condition = colData(dds)$condition
)
p <- ggplot(lib_df, aes(x=reorder(sample, counts), y=counts/1e6, fill=condition)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=cond_cols) +
  coord_flip() +
  labs(title="Library sizes — GSE201377",
       x="Sample", y="Mapped reads (millions)") +
  theme_bw(base_size=11)
ggsave(file.path(outdir, "library_sizes.pdf"), p, width=7, height=5)

cat("=== All GSE201377 plots complete ===\n")
cat("Output directory:", outdir, "\n")
