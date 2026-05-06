library(DESeq2)

raw <- read.table("/users/fyp/fyp5/project/GSE201377/featurecounts/counts_matrix.txt",
                  header=TRUE, skip=1, sep="\t", row.names=1, check.names=FALSE)
count_mat <- as.matrix(raw[, 6:ncol(raw)])
colnames(count_mat) <- gsub(".*/GSM", "GSM", colnames(count_mat))
colnames(count_mat) <- gsub("/Aligned.sortedByCoord.out.bam", "", colnames(count_mat))

condition <- ifelse(grepl("Pta", colnames(count_mat)), "Pta", "Psy")
coldata <- data.frame(condition=factor(condition, levels=c("Pta","Psy")),
                      row.names=colnames(count_mat))

dds <- DESeqDataSetFromMatrix(countData=count_mat, colData=coldata, design=~condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after filtering:", nrow(dds), "\n")

dds <- DESeq(dds)
cat("Size factors:\n"); print(sizeFactors(dds))

saveRDS(dds, "/users/fyp/fyp5/project/GSE201377/deseq2_steps/dds.rds")
cat("dds saved.\n")
