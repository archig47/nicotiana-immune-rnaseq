library(dplyr)

gff <- read.table("/users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3",
                  sep="\t", quote="", comment.char="#",
                  col.names=c("chr","source","type","start","end","score","strand","phase","attrs"))
gff <- gff[gff$type == "gene", ]
gff$gene_id   <- gsub(".*ID=([^;]+).*", "\\1", gff$attrs)
gff$gene_name <- ifelse(grepl("gene_name=", gff$attrs), gsub(".*gene_name=([^;]+).*", "\\1", gff$attrs), NA)
gff$note      <- ifelse(grepl("Note=", gff$attrs), gsub(".*Note=([^;]+).*", "\\1", gff$attrs), NA)

degs <- read.csv("/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_sig.csv")
degs_ann <- left_join(degs, gff[, c("gene_id","gene_name","note","chr","start","end","strand")], by="gene_id")

cat(sprintf("DEGs with gene_name: %d / %d\n", sum(!is.na(degs_ann$gene_name)), nrow(degs_ann)))
cat(sprintf("DEGs unannotated: %d\n", sum(is.na(degs_ann$gene_name))))

cat("\nTop 30 upregulated in Psy:\n")
up <- degs_ann[order(-degs_ann$log2FoldChange), c("gene_id","gene_name","note","log2FoldChange","padj","baseMean")]
print(head(up, 30), row.names=FALSE)

cat("\nTop 10 downregulated in Psy:\n")
down <- degs_ann[order(degs_ann$log2FoldChange), c("gene_id","gene_name","note","log2FoldChange","padj","baseMean")]
print(head(down, 10), row.names=FALSE)

write.csv(degs_ann, "/users/fyp/fyp5/project/GSE201377/deseq2_softmasked/DEGs_annotated.csv", row.names=FALSE)
cat("\nSaved: DEGs_annotated.csv\n")
