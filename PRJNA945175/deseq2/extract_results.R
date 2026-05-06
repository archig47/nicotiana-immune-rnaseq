library(DESeq2)
library(dplyr)

outdir     <- "/users/fyp/fyp5/project/PRJNA945175/deseq2"
nlrtracker <- "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.NLRtracker_classified.tsv"

cat("Loading dds...\n")
dds <- readRDS(file.path(outdir, "dds.rds"))

for (cond in c("HCRV","TSWV","TH")) {
  cat("Results:", cond, "vs CK...\n")
  res        <- results(dds, contrast=c("condition", cond, "CK"), alpha=0.05)
  res_shrunk <- lfcShrink(dds, contrast=c("condition", cond, "CK"), type="normal", res=res)
  res_df     <- as.data.frame(res_shrunk)
  res_df$gene_id <- rownames(res_df)
  write.csv(res_df, file.path(outdir, paste0("all_results_", cond, "_vs_CK.csv")), row.names=FALSE)
  degs <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  write.csv(degs, file.path(outdir, paste0("DEGs_", cond, "_vs_CK.csv")), row.names=FALSE)
  cat("  DEGs:", nrow(degs), "\n")
}

cat("Computing rlog...\n")
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file.path(outdir, "rld.rds"))

cat("Filtering NLR DEGs...\n")
nlr <- read.table(nlrtracker, header=TRUE, sep="\t")
nlr$gene_id <- gsub("-mRNA$", "", nlr$protein_id)
for (cond in c("HCRV","TSWV","TH")) {
  degs     <- read.csv(file.path(outdir, paste0("DEGs_", cond, "_vs_CK.csv")))
  nlr_degs <- degs %>% filter(gene_id %in% nlr$gene_id)
  write.csv(nlr_degs, file.path(outdir, paste0("DEGs_NLR_", cond, "_vs_CK.csv")), row.names=FALSE)
  cat("  NLR DEGs", cond, ":", nrow(nlr_degs), "\n")
}

cat("=== Results extraction complete ===\n")
