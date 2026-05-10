library(gprofiler2)

# Load PANNZER2 annotations
go_annot <- read.table("/users/fyp/fyp5/project/genome/PRR_GO_filtered.tsv",
                       header=TRUE, sep="\t", quote="")
go_annot$gene_id_clean <- sub("-mRNA$", "", go_annot$gene_id)

# Build custom GMT file
go_list <- split(go_annot$gene_id_clean, go_annot$go_term)
go_list <- go_list[sapply(go_list, length) >= 3]
cat("GO terms in custom database:", length(go_list), "\n")

# Write GMT file
gmt_path <- "/users/fyp/fyp5/project/genome/PRR_custom.gmt"
con <- file(gmt_path, "w")
for (term in names(go_list)) {
  writeLines(paste(c(term, term, go_list[[term]]), collapse="\t"), con)
}
close(con)
cat("GMT file written to", gmt_path, "\n")

# Upload to gprofiler
custom_gmt <- upload_GMT_file(gmtfile=gmt_path)
cat("Custom GMT token:", custom_gmt, "\n")

# Load DEG subsets
hcrv <- read.csv("/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs_PRR_HCRV_vs_CK.csv")
tswv <- read.csv("/users/fyp/fyp5/project/PRJNA945175/deseq2/DEGs_PRR_TSWV_vs_CK.csv")

hcrv_up <- hcrv$gene_id[hcrv$log2FoldChange > 0]
tswv_down <- tswv$gene_id[tswv$log2FoldChange < 0]
hcrv_down <- hcrv$gene_id[hcrv$log2FoldChange < 0]
tswv_specific_down <- tswv_down[!tswv_down %in% hcrv_down]
bg_genes <- unique(go_annot$gene_id_clean)

cat("HCRV_up:", length(hcrv_up), "\n")
cat("TSWV_specific_down:", length(tswv_specific_down), "\n")

# Run enrichment
run_enrich <- function(genes, label, token, bg) {
  cat("\nRunning:", label, "\n")
  res <- tryCatch(
    gost(query=genes, organism=token, custom_bg=bg,
         correction_method="fdr", significant=TRUE),
    error=function(e) { cat("Error:", e$message, "\n"); NULL }
  )
  if (is.null(res) || is.null(res$result)) {
    cat("No significant terms for", label, "\n")
    return(NULL)
  }
  df <- res$result[order(res$result$p_value),]
  cat("Significant terms:", nrow(df), "\n")
  print(head(df[, c("term_id","term_name","p_value","intersection_size","term_size")], 15))
  df
}

res_hcrv_up   <- run_enrich(hcrv_up, "HCRV_up", custom_gmt, bg_genes)
res_tswv_spec <- run_enrich(tswv_specific_down, "TSWV_specific_down", custom_gmt, bg_genes)

outdir <- "/users/fyp/fyp5/project/PRJNA945175/go_enrichment/"
dir.create(outdir, showWarnings=FALSE)
if (!is.null(res_hcrv_up))   write.csv(res_hcrv_up,   paste0(outdir, "gprofiler_HCRV_up.csv"),            row.names=FALSE)
if (!is.null(res_tswv_spec)) write.csv(res_tswv_spec, paste0(outdir, "gprofiler_TSWV_specific_down.csv"), row.names=FALSE)
cat("\nDone.\n")
# Fix save for tswv_spec
if (!is.null(res_tswv_spec)) {
  df <- res_tswv_spec[, sapply(res_tswv_spec, function(x) !is.list(x))]
  write.csv(df, paste0(outdir, "gprofiler_TSWV_specific_down.csv"), row.names=FALSE)
}
