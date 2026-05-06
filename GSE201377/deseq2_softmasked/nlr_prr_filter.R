library(dplyr)

degs <- read.csv("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_sig.csv")
all_res <- read.csv("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/all_results.csv")

nlr <- read.delim("/users/fyp/fyp5/project/genome/jiorgos/NbT2T.NLRtracker_classified.tsv", sep="\t", header=TRUE)
prr <- read.delim("/users/fyp/fyp5/project/genome/jiorgos/NbT2T.PRRtracker_classified.tsv", sep="\t", header=TRUE)

nlr$gene_id <- gsub("-mRNA$", "", nlr$protein_id)
prr$gene_id <- gsub("-mRNA$", "", prr$protein_id)

degs_nlr <- inner_join(degs, nlr, by="gene_id")
degs_prr <- inner_join(degs, prr, by="gene_id")

n_expressed <- 40944
n_degs <- nrow(degs)

fisher_test <- function(tracker_df, tracker_all, label) {
  n_tracker_expressed <- sum(tracker_all$gene_id %in% all_res$gene_id)
  mat <- matrix(c(
    nrow(tracker_df),
    n_degs - nrow(tracker_df),
    n_tracker_expressed - nrow(tracker_df),
    n_expressed - n_degs - (n_tracker_expressed - nrow(tracker_df))
  ), nrow=2)
  ft <- fisher.test(mat, alternative="greater")
  cat(sprintf("\n── %s Enrichment ──\n", label))
  cat(sprintf("  DEGs overlapping %s: %d / %d expressed\n", label, nrow(tracker_df), n_tracker_expressed))
  cat(sprintf("  Odds ratio: %.2f, p-value: %.4f\n", ft$estimate, ft$p.value))
}

cat("── NLR DEGs ──\n")
cat(sprintf("Total DEGs: %d\n", n_degs))
cat(sprintf("NLR DEGs: %d\n", nrow(degs_nlr)))
cat(sprintf("  Up in Psy: %d\n", sum(degs_nlr$log2FoldChange > 0)))
cat(sprintf("  Down in Psy: %d\n", sum(degs_nlr$log2FoldChange < 0)))
print(table(degs_nlr$type))

cat("\n── PRR DEGs ──\n")
cat(sprintf("PRR DEGs: %d\n", nrow(degs_prr)))
cat(sprintf("  Up in Psy: %d\n", sum(degs_prr$log2FoldChange > 0)))
cat(sprintf("  Down in Psy: %d\n", sum(degs_prr$log2FoldChange < 0)))
print(table(degs_prr$type))

fisher_test(degs_nlr, nlr, "NLR")
fisher_test(degs_prr, prr, "PRR")

write.csv(degs_nlr, "/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_NLR.csv", row.names=FALSE)
write.csv(degs_prr, "/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_PRR.csv", row.names=FALSE)
cat("\nSaved: DEGs_NLR.csv and DEGs_PRR.csv\n")
