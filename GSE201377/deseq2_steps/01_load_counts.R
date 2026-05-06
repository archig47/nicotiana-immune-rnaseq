library(DESeq2)

raw <- read.table("/users/fyp/fyp5/project/GSE201377/featurecounts/counts_matrix.txt",
                  header=TRUE, skip=1, sep="\t", row.names=1, check.names=FALSE)

count_mat <- as.matrix(raw[, 6:ncol(raw)])
colnames(count_mat) <- gsub(".*/GSM", "GSM", colnames(count_mat))
colnames(count_mat) <- gsub("/Aligned.sortedByCoord.out.bam", "", colnames(count_mat))

cat("Samples:\n"); print(colnames(count_mat))
cat("\nDimensions:", nrow(count_mat), "genes x", ncol(count_mat), "samples\n")
