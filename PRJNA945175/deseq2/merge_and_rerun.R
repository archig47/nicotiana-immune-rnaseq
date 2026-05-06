library(DESeq2)
library(dplyr)

outdir     <- "/users/fyp/fyp5/project/PRJNA945175/deseq2"
nlrtracker <- "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.NLRtracker_classified.tsv"
prrtracker <- "/users/fyp/fyp5/project/genome/jiorgos/NbT2T.PRRtracker_classified.tsv"

# ── Load existing 31-sample counts ────────────────────────────────────────
cat("Loading 31-sample counts...\n")
raw31 <- read.table(
  "/users/fyp/fyp5/project/PRJNA945175/featurecounts/counts_v12_multimapper.txt",
  header=TRUE, sep="\t", comment.char="#", check.names=FALSE)
counts31 <- raw31[, 7:ncol(raw31)]
rownames(counts31) <- gsub("-mRNA$", "", raw31$Geneid)
colnames(counts31) <- gsub(".*/aligned/(SRR[0-9]+)/.*", "\\1", colnames(counts31))

# ── Load SRR23875081 counts ───────────────────────────────────────────────
cat("Loading SRR23875081 counts...\n")
raw081 <- read.table(
  "/users/fyp/fyp5/project/PRJNA945175/featurecounts/counts_081_multimapper.txt",
  header=TRUE, sep="\t", comment.char="#", check.names=FALSE)
counts081 <- raw081[, 7, drop=FALSE]
rownames(counts081) <- gsub("-mRNA$", "", raw081$Geneid)
colnames(counts081) <- "SRR23875081"

# ── Merge ─────────────────────────────────────────────────────────────────
cat("Merging counts matrices...\n")
counts <- cbind(counts31, counts081)
counts <- round(counts)
cat("Total samples:", ncol(counts), "\n")
cat("Total genes:", nrow(counts), "\n")

# ── Coldata — all 32 samples ──────────────────────────────────────────────
coldata <- data.frame(
  sample    = c("SRR23875051","SRR23875052","SRR23875053","SRR23875054",
                "SRR23875055","SRR23875056","SRR23875057","SRR23875058",
                "SRR23875060","SRR23875061","SRR23875062","SRR23875063",
                "SRR23875064","SRR23875065","SRR23875066","SRR23875067",
                "SRR23875068","SRR23875069","SRR23875070","SRR23875072",
                "SRR23875073","SRR23875074","SRR23875075","SRR23875077",
                "SRR23875079","SRR23875080","SRR23875081","SRR23875082",
                "SRR23875083","SRR23875084","SRR23875085","SRR23875086"),
  condition = c("HCRV","TSWV","TSWV","TSWV","TSWV","TSWV","CK","CK",
                "CK","CK","CK","TSWV","CK","CK","CK","TH","TH","TH",
                "TH","TH","TH","TSWV","TH","HCRV","HCRV","HCRV","HCRV",
                "HCRV","HCRV","HCRV","TSWV","TSWV"),
  timepoint = c("1","14","14","14","7","7","14","14",
                "7","7","7","7","1","1","1","14","14","14",
                "7","7","1","1","1","14","14","7","7",
                "7","1","1","1","1"),
  stringsAsFactors = FALSE
)
coldata$condition <- factor(coldata$condition, levels=c("CK","HCRV","TSWV","TH"))
coldata$timepoint <- factor(coldata$timepoint, levels=c("1","7","14"))
rownames(coldata) <- coldata$sample
counts <- counts[, coldata$sample]

# ── DESeq2 ────────────────────────────────────────────────────────────────
cat("Running DESeq2 with all 32 samples...\n")
dds  <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~timepoint+condition)
keep <- rowSums(counts(dds) >= 10) >= 2
dds  <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")
dds  <- DESeq(dds)
saveRDS(dds, file.path(outdir, "dds_32.rds"))
cat("dds_32 saved\n")

# ── Results ───────────────────────────────────────────────────────────────
for (cond in c("HCRV","TSWV","TH")) {
  cat("Results:", cond, "vs CK...\n")
  res        <- results(dds, contrast=c("condition", cond, "CK"), alpha=0.05)
  res_shrunk <- lfcShrink(dds, contrast=c("condition", cond, "CK"), type="normal", res=res)
  res_df     <- as.data.frame(res_shrunk)
  res_df$gene_id <- rownames(res_df)
  write.csv(res_df, file.path(outdir, paste0("all_results32_", cond, "_vs_CK.csv")), row.names=FALSE)
  degs <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
  write.csv(degs, file.path(outdir, paste0("DEGs32_", cond, "_vs_CK.csv")), row.names=FALSE)
  cat("  DEGs:", nrow(degs), "\n")
}

# ── rlog ──────────────────────────────────────────────────────────────────
cat("Computing rlog...\n")
rld <- rlog(dds, blind=FALSE)
saveRDS(rld, file.path(outdir, "rld_32.rds"))

# ── NLR + PRR overlap ─────────────────────────────────────────────────────
cat("Filtering NLR and PRR DEGs...\n")
nlr <- read.table(nlrtracker, header=TRUE, sep="\t")
nlr$gene_id <- gsub("-mRNA$", "", nlr$protein_id)
prr <- read.table(prrtracker, header=TRUE, sep="\t")
prr$gene_id <- gsub("-mRNA$", "", prr$protein_id)

for (cond in c("HCRV","TSWV","TH")) {
  degs <- read.csv(file.path(outdir, paste0("DEGs32_", cond, "_vs_CK.csv")))
  nlr_degs <- degs %>% filter(gene_id %in% nlr$gene_id)
  prr_degs <- degs %>% filter(gene_id %in% prr$gene_id)
  write.csv(nlr_degs, file.path(outdir, paste0("DEGs32_NLR_", cond, "_vs_CK.csv")), row.names=FALSE)
  write.csv(prr_degs, file.path(outdir, paste0("DEGs32_PRR_", cond, "_vs_CK.csv")), row.names=FALSE)
  cat("  NLR DEGs", cond, ":", nrow(nlr_degs), "\n")
  cat("  PRR DEGs", cond, ":", nrow(prr_degs), "\n")
}

cat("=== Complete ===\n")
