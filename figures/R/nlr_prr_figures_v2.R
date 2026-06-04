#!/usr/bin/env Rscript
# ============================================================
# NLR and PRR DEG visualisations — refined subclass from InterPro
# 1. Heatmaps ordered by subclass
# 2. Abundance vs LFC scatter plots by timepoint
# ============================================================

suppressPackageStartupMessages({
  library(pheatmap)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

# ── Paths ──────────────────────────────────────────────────
EXPR_MAT  <- "/users/fyp/fyp5/project/PRJNA945175/clust/expression_matrix_tp.tsv"
LANDSCAPE <- "/users/fyp/fyp5/project/PRJNA945175/deseq2/nlr_prr_full_landscape.xlsx"
PRR_CLUST <- "/users/fyp/fyp5/project/PRJNA945175/clust/clust_output_prr_degs/Clusters_Objects.tsv"
DEG_DIR   <- "/users/fyp/fyp5/project/PRJNA945175/deseq2"
OUTDIR    <- "/users/fyp/fyp5/project/PRJNA945175/figures"
dir.create(OUTDIR, showWarnings=FALSE)

# ── Load expression matrix ─────────────────────────────────
cat("Loading expression matrix...\n")
expr <- read.table(EXPR_MAT, header=TRUE, sep="\t", row.names=1, check.names=FALSE)
col_order <- c("CK_1","CK_7","CK_14",
               "HCRV_1","HCRV_7","HCRV_14",
               "TH_1","TH_7","TH_14",
               "TSWV_1","TSWV_7","TSWV_14")
expr <- expr[, col_order]

# ── Load landscape ─────────────────────────────────────────
nlr <- read_xlsx(LANDSCAPE, sheet="NLR_landscape")
prr <- read_xlsx(LANDSCAPE, sheet="PRR_landscape")
nlr_degs <- nlr %>% filter(is_DEG == TRUE)
prr_degs <- prr %>% filter(is_DEG == TRUE)

# ── NLR subclass — refined using domain_arch_simple + InterPro ──
nlr_degs <- nlr_degs %>%
  mutate(subclass = case_when(
    # ROQ1-type TNL — biologically distinct, worth separating
    grepl("WHD_ROQ1", interpro_names)                          ~ "TNL-ROQ1",
    # TIR-containing = TNL
    grepl("T.*NL|^TN|TNT", domain_arch_simple)                 ~ "TNL",
    # CC-containing with CJID = CNL-CJID
    grepl("C-JID", interpro_names)                             ~ "CNL-CJID",
    # RPW8 domain
    grepl("Powdery_mildew", interpro_names) | type == "RPW8"   ~ "RPW8",
    # Jacalin-integrated NLR
    grepl("Jacalin", interpro_names)                           ~ "NLR-ID (Jacalin)",
    # CC-containing = CNL (includes Rx-type)
    grepl("C.*NL|^CN|^NC|BCN|BNL|CNCL|CNOL|ESF2", domain_arch_simple) |
    grepl("RX-like_CC|Rx_N|ESF2", interpro_names)             ~ "CNL",
    # NL = NB-ARC+LRR, no N-terminal confirmed
    grepl("^NL$|^NLJ$|^CNOL$", domain_arch_simple)            ~ "NLR (no N-term)",
    # N only = NB-ARC only
    grepl("^N$", domain_arch_simple)                           ~ "NLR (NB-ARC only)",
    # CCX, TX classes from tracker
    type == "CCX"                                              ~ "CCX",
    type == "TX"                                               ~ "TX",
    type == "degenerate_NLR"                                   ~ "Degenerate",
    TRUE                                                       ~ "Other"
  ))

cat("NLR subclass breakdown:\n")
print(sort(table(nlr_degs$subclass), decreasing=TRUE))

# ── PRR subclass — refined using matched_rules + InterPro ──
prr_degs <- prr_degs %>%
  mutate(subclass = case_when(
    # LysM-RLK — chitin/MAMP perception
    grepl("LysM", interpro_names)                              ~ "LysM-RLK",
    # B-lectin RLK
    grepl("BLEC", matched_rules) |
    grepl("Bulb-type_lectin|ConA-like", interpro_names)        ~ "LecRK",
    # WAK — wall-associated kinase
    grepl("WAK", matched_rules) |
    grepl("WAK_GUB|WAK_assoc", interpro_names)                 ~ "WAK",
    # Malectin-RLK
    grepl("Malectin", matched_rules) |
    grepl("Malectin", interpro_names)                          ~ "Malectin-RLK",
    # BRI1-type LRR-RLK
    grepl("BRI1_island", interpro_names)                       ~ "BRI1-type LRR-RLK",
    # GNK2-type
    grepl("GNK2|Gnk2", interpro_names)                        ~ "GNK2-RLK",
    # LRR-RLK (generic LRR + kinase)
    grepl("LRR", matched_rules) & type == "RLK"               ~ "LRR-RLK",
    # RLK-other (kinase but no identified ectodomain)
    type == "RLK"                                              ~ "RLK-other",
    # RLP
    type == "RLP"                                              ~ "RLP",
    # Secreted
    type == "Secreted"                                         ~ "Secreted",
    TRUE                                                       ~ "Other"
  ))

cat("\nPRR subclass breakdown:\n")
print(sort(table(prr_degs$subclass), decreasing=TRUE))

# ── Load PRR clusters ──────────────────────────────────────
clust_raw <- read.table(PRR_CLUST, header=FALSE, sep="\t",
                        fill=TRUE, stringsAsFactors=FALSE)
clust_names <- as.character(clust_raw[1, ])
clust_map <- data.frame(gene_id=character(), cluster=character(),
                        stringsAsFactors=FALSE)
for (i in seq_along(clust_names)) {
  cname <- clust_names[i]
  if (nchar(cname) > 0 && grepl("^C", cname)) {
    cname_short <- gsub(" .*", "", cname)
    genes <- as.character(clust_raw[-c(1,2), i])
    genes <- genes[nchar(genes) > 0 & grepl("NbT2T", genes)]
    if (length(genes) > 0)
      clust_map <- rbind(clust_map,
                         data.frame(gene_id=genes, cluster=cname_short))
  }
}
prr_degs <- prr_degs %>%
  left_join(clust_map, by="gene_id") %>%
  mutate(cluster = ifelse(is.na(cluster), "Unclustered", cluster))

# ── Load per-contrast DESeq2 results ──────────────────────
cat("\nLoading DESeq2 results...\n")
contrasts <- c("HCRV_1","HCRV_7","HCRV_14",
               "TH_1","TH_7","TH_14",
               "TSWV_1","TSWV_7","TSWV_14")
deg_list <- list()
for (ct in contrasts) {
  f <- file.path(DEG_DIR, paste0("results_", ct, "dpi_vs_CK.csv"))
  if (file.exists(f)) {
    d <- read.csv(f, row.names=1) %>%
      filter(!is.na(padj), padj < 0.05) %>%
      mutate(gene_id   = rownames(.),
             contrast  = ct,
             timepoint = sub(".*_(\\d+)$", "\\1 dpi", ct),
             condition = sub("_\\d+$", "", ct))
    deg_list[[ct]] <- d
  }
}
all_degs <- bind_rows(deg_list)
cat("Total per-contrast DEG rows:", nrow(all_degs), "\n")

# ── Colour palettes ────────────────────────────────────────
nlr_subclass_cols <- c(
  "CNL"              = "#2166AC",
  "CNL-CJID"         = "#6BAED6",
  "TNL"              = "#1A9641",
  "TNL-ROQ1"         = "#74C476",
  "RPW8"             = "#E66101",
  "NLR (no N-term)"  = "#9ECAE1",
  "NLR (NB-ARC only)"= "#DEEBF7",
  "NLR-ID (Jacalin)" = "#FDD0A2",
  "CCX"              = "#D01C8B",
  "TX"               = "#FDB863",
  "Degenerate"       = "#B2ABD2",
  "Other"            = "#CCCCCC"
)

prr_subclass_cols <- c(
  "LRR-RLK"       = "#1A9641",
  "LecRK"         = "#74C476",
  "WAK"           = "#FEE08B",
  "Malectin-RLK"  = "#FDAE61",
  "LysM-RLK"      = "#D73027",
  "BRI1-type LRR-RLK" = "#ABD9E9",
  "GNK2-RLK"      = "#F46D43",
  "RLK-other"     = "#D9EF8B",
  "RLP"           = "#74ADD1",
  "Secreted"      = "#C6DBEF",
  "Other"         = "#CCCCCC"
)

dir_cols   <- c("Up"="#D73027","Down"="#4575B4")
clust_cols <- c("C0"="#D7191C","C1"="#2C7BB6","Unclustered"="#DDDDDD")

heatmap_cols <- colorRampPalette(
  c("#4575B4","#91BFDB","#E0F3F8","white","#FEE090","#FC8D59","#D73027"))(100)

col_ann <- data.frame(
  Condition = c(rep("CK",3),rep("HCRV",3),rep("TH",3),rep("TSWV",3)),
  Timepoint = rep(c("1 dpi","7 dpi","14 dpi"),4),
  row.names = col_order
)
cond_cols <- c("CK"="#888888","HCRV"="#E69F00","TH"="#56B4E9","TSWV"="#D55E00")
tp_cols   <- c("1 dpi"="#FEE08B","7 dpi"="#FDAE61","14 dpi"="#F46D43")
col_ann_colors <- list(Condition=cond_cols, Timepoint=tp_cols)

# ════════════════════════════════════════════════════════════
# HEATMAP 1 — NLR DEGs
# ════════════════════════════════════════════════════════════
cat("\n-- NLR Heatmap --\n")

nlr_level_order <- c("CNL","CNL-CJID","TNL","TNL-ROQ1","RPW8",
                     "NLR (no N-term)","NLR (NB-ARC only)",
                     "NLR-ID (Jacalin)","CCX","TX","Degenerate","Other")

nlr_ordered <- nlr_degs %>%
  mutate(subclass = factor(subclass, levels=nlr_level_order)) %>%
  arrange(subclass, desc(mean_abs_LFC))

nlr_ids  <- nlr_ordered$gene_id
nlr_expr <- expr[rownames(expr) %in% nlr_ids, , drop=FALSE]
nlr_expr <- nlr_expr[nlr_ids[nlr_ids %in% rownames(nlr_expr)], ]
cat("NLR genes in matrix:", nrow(nlr_expr), "\n")

nlr_row_ann <- nlr_ordered %>%
  filter(gene_id %in% rownames(nlr_expr)) %>%
  select(gene_id, subclass, direction) %>%
  tibble::column_to_rownames("gene_id") %>%
  rename(Subclass=subclass, Direction=direction)

nlr_ann_colors <- list(
  Subclass  = nlr_subclass_cols[names(nlr_subclass_cols) %in% unique(as.character(nlr_row_ann$Subclass))],
  Direction = dir_cols
)

nlr_mat <- t(scale(t(as.matrix(nlr_expr))))
nlr_mat[nlr_mat >  3] <-  3
nlr_mat[nlr_mat < -3] <- -3

# Compute gaps between subclass groups
nlr_gap_counts <- nlr_ordered %>%
  filter(gene_id %in% rownames(nlr_expr)) %>%
  group_by(subclass) %>%
  summarise(n=n(), .groups="drop") %>%
  arrange(match(subclass, nlr_level_order))
nlr_gaps <- cumsum(nlr_gap_counts$n)[-nrow(nlr_gap_counts)]

pheatmap(nlr_mat,
  color             = heatmap_cols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_row    = nlr_row_ann,
  annotation_col    = col_ann,
  annotation_colors = c(nlr_ann_colors, col_ann_colors),
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  fontsize_col      = 9,
  fontsize          = 9,
  border_color      = NA,
  gaps_row          = nlr_gaps,
  gaps_col          = c(3, 6, 9),
  main              = paste0("NLR DEGs (n=", nrow(nlr_mat), ") ordered by subclass"),
  filename          = file.path(OUTDIR, "heatmap_nlr_degs.pdf"),
  width=8, height=10
)
cat("Saved: heatmap_nlr_degs.pdf\n")

# ════════════════════════════════════════════════════════════
# HEATMAP 2 — PRR clustered (C0+C1)
# ════════════════════════════════════════════════════════════
cat("\n-- PRR Heatmap (clustered) --\n")

prr_level_order <- c("LRR-RLK","LecRK","WAK","Malectin-RLK",
                     "LysM-RLK","BRI1-type LRR-RLK","GNK2-RLK",
                     "RLK-other","RLP","Secreted","Other")

prr_clustered <- prr_degs %>%
  filter(cluster %in% c("C0","C1")) %>%
  mutate(subclass = factor(subclass, levels=prr_level_order),
         cluster  = factor(cluster, levels=c("C0","C1"))) %>%
  arrange(cluster, subclass, desc(mean_abs_LFC))

prr_ids  <- prr_clustered$gene_id
prr_expr <- expr[rownames(expr) %in% prr_ids, , drop=FALSE]
prr_expr <- prr_expr[prr_ids[prr_ids %in% rownames(prr_expr)], ]
cat("PRR clustered genes in matrix:", nrow(prr_expr), "\n")

prr_row_ann <- prr_clustered %>%
  filter(gene_id %in% rownames(prr_expr)) %>%
  select(gene_id, cluster, subclass, direction) %>%
  tibble::column_to_rownames("gene_id") %>%
  rename(Cluster=cluster, Subclass=subclass, Direction=direction)

prr_ann_colors <- list(
  Cluster   = clust_cols[c("C0","C1")],
  Subclass  = prr_subclass_cols[names(prr_subclass_cols) %in% unique(as.character(prr_row_ann$Subclass))],
  Direction = dir_cols
)

prr_mat <- t(scale(t(as.matrix(prr_expr))))
prr_mat[prr_mat >  3] <-  3
prr_mat[prr_mat < -3] <- -3

prr_c0_n <- sum(prr_clustered$cluster[prr_clustered$gene_id %in% rownames(prr_expr)] == "C0")

pheatmap(prr_mat,
  color             = heatmap_cols,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  annotation_row    = prr_row_ann,
  annotation_col    = col_ann,
  annotation_colors = c(prr_ann_colors, col_ann_colors),
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  fontsize_col      = 9,
  fontsize          = 9,
  border_color      = NA,
  gaps_row          = prr_c0_n,
  gaps_col          = c(3, 6, 9),
  main              = paste0("PRR co-expression clusters C0+C1 (n=", nrow(prr_mat), ")"),
  filename          = file.path(OUTDIR, "heatmap_prr_clustered.pdf"),
  width=8, height=8
)
cat("Saved: heatmap_prr_clustered.pdf\n")

# ════════════════════════════════════════════════════════════
# SCATTER PLOTS — Abundance vs LFC by timepoint
# ════════════════════════════════════════════════════════════
cat("\n-- Scatter plots --\n")

make_scatter_pdf <- function(degs_landscape, all_degs_df,
                             subclass_cols, outfile, label_n=5) {

  dat <- all_degs_df %>%
    inner_join(degs_landscape %>% select(gene_id, subclass), by="gene_id") %>%
    mutate(
      timepoint      = factor(timepoint, levels=c("1 dpi","7 dpi","14 dpi")),
      log10_baseMean = log10(baseMean + 1),
      subclass       = factor(subclass, levels=names(subclass_cols))
    )

  if (nrow(dat) == 0) { cat("No data\n"); return(NULL) }

  pdf(outfile, width=14, height=5)
  for (tp in c("1 dpi","7 dpi","14 dpi")) {
    d <- dat %>% filter(timepoint == tp)
    if (nrow(d) == 0) next

    top_genes <- d %>%
      filter(baseMean > 100) %>%
      arrange(desc(abs(log2FoldChange))) %>%
      slice_head(n=label_n)

    p <- ggplot(d, aes(x=log10_baseMean, y=log2FoldChange, colour=subclass)) +
      geom_hline(yintercept=0, linetype="dashed", colour="grey60", linewidth=0.4) +
      geom_hline(yintercept=c(-1,1), linetype="dotted", colour="grey80", linewidth=0.3) +
      geom_point(alpha=0.7, size=1.8) +
      geom_text_repel(data=top_genes,
        aes(label=gene_id), size=2.2, max.overlaps=15,
        segment.size=0.3, segment.colour="grey40",
        box.padding=0.3) +
      scale_colour_manual(values=subclass_cols, name="Subclass", drop=FALSE) +
      scale_x_continuous(name=expression(log[10](baseMean + 1)),
                         limits=c(0, NA)) +
      scale_y_continuous(name=expression(log[2]~Fold~Change)) +
      facet_wrap(~condition, nrow=1, scales="free_x") +
      ggtitle(paste0("Abundance vs Fold Change — ", tp)) +
      theme_bw(base_size=9) +
      theme(
        strip.background  = element_rect(fill="grey92", colour=NA),
        strip.text        = element_text(face="bold"),
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(colour="grey95"),
        legend.position   = "right",
        legend.key.size   = unit(0.4,"cm"),
        plot.title        = element_text(face="bold", size=11),
        axis.title        = element_text(size=9)
      )
    print(p)
  }
  dev.off()
  cat("Saved:", outfile, "\n")
}

# NLR scatter
make_scatter_pdf(
  degs_landscape = nlr_degs,
  all_degs_df    = all_degs,
  subclass_cols  = nlr_subclass_cols,
  outfile        = file.path(OUTDIR, "scatter_nlr_abundance_lfc.pdf"),
  label_n        = 5
)

# PRR scatter
make_scatter_pdf(
  degs_landscape = prr_degs,
  all_degs_df    = all_degs,
  subclass_cols  = prr_subclass_cols,
  outfile        = file.path(OUTDIR, "scatter_prr_abundance_lfc.pdf"),
  label_n        = 5
)

cat("\nAll figures saved to:", OUTDIR, "\n")
