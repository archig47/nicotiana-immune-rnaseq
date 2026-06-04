suppressPackageStartupMessages({
  library(DESeq2)
  library(readr)
  library(readxl)
  library(dplyr)
  library(openxlsx)
})

CLEAN_DIR  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2_clean"
RESULTS    <- file.path(CLEAN_DIR, "results")
LANDSCAPE  <- "/users/fyp/fyp5/project/PRJNA945175/deseq2/nlr_prr_full_landscape.xlsx"
OUTPUT     <- file.path(CLEAN_DIR, "nlr_prr", "nlr_prr_clean_summary.xlsx")

PADJ_THRESH <- 0.05
LFC_THRESH  <- 1

# ── 1. Load landscape gene IDs ────────────────────────────────
cat("Loading landscape...\n")
nlr_landscape <- read_xlsx(LANDSCAPE, sheet="NLR_landscape")
prr_landscape <- read_xlsx(LANDSCAPE, sheet="PRR_landscape")

nlr_ids <- nlr_landscape$gene_id
prr_ids <- prr_landscape$gene_id
cat(sprintf("NLR genes in genome: %d\n", length(nlr_ids)))
cat(sprintf("PRR genes in genome: %d\n", length(prr_ids)))

# ── 2. Load all 9 contrast results ────────────────────────────
cat("\nLoading contrast results...\n")
contrasts <- list(
  HCRV_1  = "results_HCRV_1dpi_vs_CK.csv",
  HCRV_7  = "results_HCRV_7dpi_vs_CK.csv",
  HCRV_14 = "results_HCRV_14dpi_vs_CK.csv",
  TSWV_1  = "results_TSWV_1dpi_vs_CK.csv",
  TSWV_7  = "results_TSWV_7dpi_vs_CK.csv",
  TSWV_14 = "results_TSWV_14dpi_vs_CK.csv",
  TH_1    = "results_TH_1dpi_vs_CK.csv",
  TH_7    = "results_TH_7dpi_vs_CK.csv",
  TH_14   = "results_TH_14dpi_vs_CK.csv"
)

all_results <- lapply(names(contrasts), function(k) {
  df <- read_csv(file.path(RESULTS, contrasts[[k]]), show_col_types=FALSE)
  df$contrast <- k
  df$virus    <- sub("_.*", "", k)
  df$timepoint <- as.integer(sub(".*_(\\d+)$", "\\1", k))
  df
})
all_res <- bind_rows(all_results)

# ── 3. Identify DEGs ──────────────────────────────────────────
all_res <- all_res %>%
  mutate(
    is_DEG = !is.na(padj) & padj < PADJ_THRESH & abs(log2FoldChange) > LFC_THRESH,
    direction = case_when(
      is_DEG & log2FoldChange >  LFC_THRESH ~ "up",
      is_DEG & log2FoldChange < -LFC_THRESH ~ "down",
      TRUE ~ NA_character_
    )
  )

# ── 4. Load normalised counts for expression check ───────────
cat("Loading normalised counts...\n")
dds        <- readRDS(file.path(CLEAN_DIR, "dds_clean.rds"))
norm_counts <- counts(dds, normalized=TRUE)

# A gene is expressed if norm counts >= 10 in >= 2 samples across ANY group
cat("Computing expression status...\n")

gene_ids_to_check <- unique(c(nlr_ids, prr_ids))
gene_ids_to_check <- gene_ids_to_check[gene_ids_to_check %in% rownames(norm_counts)]

expressed_any <- setNames(
  apply(norm_counts[gene_ids_to_check, ], 1, function(x) {
    # across all samples regardless of group
    sum(x >= 10) >= 2
  }),
  gene_ids_to_check
)

# ── 5. Build NLR summary table ────────────────────────────────
cat("Building NLR summary...\n")

build_summary <- function(gene_ids, landscape_df, class_col) {
  # For each gene: which contrasts is it a DEG in?
  gene_res <- all_res %>%
    filter(gene_id %in% gene_ids) %>%
    select(gene_id, contrast, virus, timepoint, log2FoldChange, padj, is_DEG, direction)

  # Per-gene summary
  per_gene <- gene_res %>%
    group_by(gene_id) %>%
    summarise(
      n_contrasts_DEG  = sum(is_DEG, na.rm=TRUE),
      is_any_DEG       = any(is_DEG, na.rm=TRUE),
      contrasts_DEG    = paste(contrast[is_DEG], collapse="; "),
      mean_LFC_DEG     = ifelse(any(is_DEG, na.rm=TRUE),
                                round(mean(log2FoldChange[is_DEG], na.rm=TRUE), 3),
                                NA_real_),
      max_abs_LFC      = round(max(abs(log2FoldChange), na.rm=TRUE), 3),
      overall_direction = ifelse(any(is_DEG, na.rm=TRUE),
                                ifelse(mean(log2FoldChange[is_DEG], na.rm=TRUE) > 0,
                                       "up", "down"),
                                NA_character_),
      peak_contrast    = ifelse(any(is_DEG, na.rm=TRUE),
                                contrast[is_DEG][which.max(abs(log2FoldChange[is_DEG]))],
                                NA_character_),
      .groups="drop"
    ) %>%
    rename(direction = overall_direction)

  # Add expressed status — gene expressed if >= 10 counts in >= 2 samples across any group
  per_gene$is_expressed <- per_gene$gene_id %in% names(expressed_any)[expressed_any]

  # Merge with landscape — rename direction before join to avoid conflict
  per_gene <- per_gene %>% rename(deg_direction = direction)

  result <- landscape_df %>%
    select(-any_of(c("direction", "is_DEG", "mean_LFC", "n_sig",
                     "mean_baseMean", "mean_abs_LFC", "peak_contrast"))) %>%
    left_join(per_gene, by="gene_id") %>%
    mutate(
      is_any_DEG      = ifelse(is.na(is_any_DEG), FALSE, is_any_DEG),
      is_expressed    = ifelse(is.na(is_expressed), FALSE, is_expressed),
      n_contrasts_DEG = ifelse(is.na(n_contrasts_DEG), 0L, n_contrasts_DEG)
    ) %>%
    arrange(desc(is_any_DEG), desc(n_contrasts_DEG))

  result
}

nlr_summary <- build_summary(nlr_ids, nlr_landscape, "nlr_class")
prr_summary <- build_summary(prr_ids, prr_landscape, "prr_class")

cat(sprintf("NLR DEGs: %d / %d\n", sum(nlr_summary$is_any_DEG), nrow(nlr_summary)))
cat(sprintf("PRR DEGs: %d / %d\n", sum(prr_summary$is_any_DEG), nrow(prr_summary)))

# ── 6. Per-contrast NLR/PRR DEG counts ───────────────────────
nlr_per_contrast <- all_res %>%
  filter(gene_id %in% nlr_ids, is_DEG) %>%
  group_by(contrast, virus, timepoint, direction) %>%
  summarise(n=n(), .groups="drop") %>%
  tidyr::pivot_wider(names_from=direction, values_from=n, values_fill=0) %>%
  mutate(total = up + down) %>%
  arrange(virus, timepoint)

prr_per_contrast <- all_res %>%
  filter(gene_id %in% prr_ids, is_DEG) %>%
  group_by(contrast, virus, timepoint, direction) %>%
  summarise(n=n(), .groups="drop") %>%
  tidyr::pivot_wider(names_from=direction, values_from=n, values_fill=0) %>%
  mutate(total = up + down) %>%
  arrange(virus, timepoint)

# ── 7. Write Excel ────────────────────────────────────────────
cat("Writing Excel...\n")
wb <- createWorkbook()

# Styles
header_style <- createStyle(
  fontColour="#FFFFFF", bgFill="#1F3864",
  textDecoration="bold", halign="center", fontSize=10
)
deg_style_up   <- createStyle(bgFill="#FCE4D6")
deg_style_down <- createStyle(bgFill="#DDEEFF")
deg_style      <- createStyle(bgFill="#E2EFDA")
normal_style   <- createStyle(fontSize=9)

write_sheet <- function(wb, sheet_name, df, deg_col="is_any_DEG", dir_col="deg_direction") {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  addStyle(wb, sheet_name, header_style,
           rows=1, cols=1:ncol(df), gridExpand=TRUE)
  setColWidths(wb, sheet_name, cols=1:ncol(df), widths="auto")
  freezePane(wb, sheet_name, firstRow=TRUE)

  if (!is.null(deg_col) && deg_col %in% names(df)) {
    for (i in 1:nrow(df)) {
      if (!is.na(df[[deg_col]][i]) && df[[deg_col]][i] == TRUE) {
        dir <- if (!is.null(dir_col) && dir_col %in% names(df)) df[[dir_col]][i] else NA
        s <- if (!is.na(dir) && dir == "up") deg_style_up else
             if (!is.na(dir) && dir == "down") deg_style_down else deg_style
        addStyle(wb, sheet_name, s, rows=i+1, cols=1:ncol(df), gridExpand=TRUE)
      }
    }
  }
}

write_sheet(wb, "NLR_all",          nlr_summary)
write_sheet(wb, "PRR_all",          prr_summary)
write_sheet(wb, "NLR_DEGs_only",    filter(nlr_summary, is_any_DEG==TRUE))
write_sheet(wb, "PRR_DEGs_only",    filter(prr_summary, is_any_DEG==TRUE))
write_sheet(wb, "NLR_per_contrast", nlr_per_contrast, deg_col=NULL)
write_sheet(wb, "PRR_per_contrast", prr_per_contrast, deg_col=NULL)

# Summary sheet
summary_df <- data.frame(
  Category = c(
    "Total NLR genes", "NLR expressed (any condition)", "NLR DEGs (any contrast)",
    "NLR DEGs upregulated", "NLR DEGs downregulated",
    "Total PRR genes", "PRR expressed (any condition)", "PRR DEGs (any contrast)",
    "PRR DEGs upregulated", "PRR DEGs downregulated"
  ),
  Count = c(
    nrow(nlr_summary),
    sum(nlr_summary$is_expressed),
    sum(nlr_summary$is_any_DEG),
    sum(nlr_summary$is_any_DEG & nlr_summary$deg_direction == "up",   na.rm=TRUE),
    sum(nlr_summary$is_any_DEG & nlr_summary$deg_direction == "down", na.rm=TRUE),
    nrow(prr_summary),
    sum(prr_summary$is_expressed),
    sum(prr_summary$is_any_DEG),
    sum(prr_summary$is_any_DEG & prr_summary$deg_direction == "up",   na.rm=TRUE),
    sum(prr_summary$is_any_DEG & prr_summary$deg_direction == "down", na.rm=TRUE)
  )
)
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)
addStyle(wb, "Summary", header_style, rows=1, cols=1:2, gridExpand=TRUE)
setColWidths(wb, "Summary", cols=1:2, widths=c(40, 15))

saveWorkbook(wb, OUTPUT, overwrite=TRUE)
cat(sprintf("\nSaved: %s\n", OUTPUT))
cat("\n=== Summary ===\n")
print(summary_df)
cat("\nNLR per contrast:\n")
print(nlr_per_contrast)
cat("\nPRR per contrast:\n")
print(prr_per_contrast)
