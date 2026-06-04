# =============================================================================
# select_marker_genes.R
# =============================================================================
# Purpose: Identify high-confidence viral-responsive marker gene candidates
#          from DESeq2 timepoint-specific contrasts and clust expression
#          clusters, using a combined abundance-activation scoring approach.
#
# Author:  Archita Gupta, Imperial College London, Kourelis Lab
# Date:    May 2025
#
# Usage:
#   Rscript select_marker_genes.R
#
# Input files required:
#   - high_conf_core_v2.csv        : High-confidence gene set from clust
#                                    abundance filtering (produced by
#                                    abundance_scatter_v2.R)
#   - NbT2T.gene_names.csv         : Gene annotations including arabidopsis
#                                    orthologs, NLR and PRR classifications
#
# Output files:
#   - marker_genes_score3.csv      : Full annotated table of selected genes
#   - score3_up_genes.txt          : Gene ID list for promoter extraction
#   - score_distribution_table.csv : Score cutoff analysis table for methods
#
# Dependencies: dplyr, readr, writexl
# =============================================================================

library(dplyr)
library(readr)
library(writexl)

# =============================================================================
# PARAMETERS — modify these to change analysis behaviour
# =============================================================================

# Score threshold for gene selection
# Score = log10(mean_baseMean + 1) x mean_abs_LFC
# Higher = more stringent. 3.0 is the natural inflection point in this dataset.
# Change to 3.5 for more stringent, 2.5 for more permissive.
SCORE_THRESHOLD <- 3.0

# Direction filter: "Up", "Down", or "both"
# "Up" selects genes upregulated relative to mock (CK)
# Change to "Down" to analyse suppressed genes, "both" for all
DIRECTION <- "Up"

# Paths — update if running on a different system
CLUST_DIR  <- "/users/fyp/fyp5/project/PRJNA945175/clust"
GENOME_DIR <- "/users/fyp/fyp5/project/genome"
OUT_DIR    <- "/users/fyp/fyp5/project/PRJNA945175/promoters"

# =============================================================================
# LOAD DATA
# =============================================================================
cat("Loading high-confidence core genes...\n")
hc <- read_csv(file.path(CLUST_DIR, "high_conf_core_v2.csv"),
               show_col_types = FALSE)
cat(sprintf("  Total high-confidence genes: %d\n", nrow(hc)))

cat("Loading gene annotations...\n")
gene_names <- read_csv(file.path(GENOME_DIR, "NbT2T.gene_names.csv"),
                       show_col_types = FALSE) %>%
  select(gene_id, name, description, arabidopsis_names, nlr_class, prr_class)

# =============================================================================
# SCORE CALCULATION
# Score = log10(mean_baseMean + 1) x mean_abs_LFC
#
# Rationale:
#   - log10(mean_baseMean + 1): rewards abundance on a log scale so that
#     differences between lowly expressed genes (100 vs 200 counts) are
#     weighted less than differences at high expression (1000 vs 2000 counts)
#   - mean_abs_LFC: rewards consistent fold change across conditions
#   - Multiplying both penalises genes that excel on only one axis —
#     a gene must be both abundant AND strongly changed to score highly
# =============================================================================
cat("\nCalculating abundance-activation scores...\n")

all_up <- hc %>%
  filter(direction == DIRECTION) %>%
  mutate(
    # Explicit score recalculation from components for transparency
    score_calc = log10(mean_baseMean + 1) * mean_abs_LFC
  ) %>%
  arrange(desc(score_calc))

cat(sprintf("  %s-regulated genes: %d\n", DIRECTION, nrow(all_up)))
cat(sprintf("  Score range: %.3f to %.3f\n",
            min(all_up$score_calc), max(all_up$score_calc)))
cat(sprintf("  Score mean: %.3f | median: %.3f\n",
            mean(all_up$score_calc), median(all_up$score_calc)))

# =============================================================================
# SCORE DISTRIBUTION TABLE
# Shows gene counts and gene properties at successive score cutoffs
# Use this table to justify your chosen threshold in your methods section
# =============================================================================
cat("\nScore distribution analysis:\n")
cat(sprintf("  Score > 5 : %d genes\n", sum(all_up$score_calc > 5)))
cat(sprintf("  Score 4-5 : %d genes\n",
            sum(all_up$score_calc >= 4 & all_up$score_calc < 5)))
cat(sprintf("  Score 3-4 : %d genes\n",
            sum(all_up$score_calc >= 3 & all_up$score_calc < 4)))
cat(sprintf("  Score 2-3 : %d genes\n",
            sum(all_up$score_calc >= 2 & all_up$score_calc < 3)))

# Check rank at key score thresholds
for (thresh in c(5, 4, 3.5, 3, 2.5)) {
  n <- sum(all_up$score_calc > thresh)
  if (n > 0 & n <= nrow(all_up)) {
    boundary <- all_up %>% filter(score_calc > thresh) %>% slice_tail(n=1)
    cat(sprintf("  Score > %.1f: %d genes | boundary gene baseMean=%.0f, LFC=%.2f\n",
                thresh, n,
                boundary$mean_baseMean,
                boundary$mean_abs_LFC))
  }
}

# Save distribution table for methods writeup
dist_table <- data.frame(
  score_threshold = c(5, 4, 3.5, 3, 2.5),
  n_genes = sapply(c(5, 4, 3.5, 3, 2.5),
                   function(t) sum(all_up$score_calc > t)),
  median_baseMean = sapply(c(5, 4, 3.5, 3, 2.5),
                           function(t) round(median(
                             all_up$mean_baseMean[all_up$score_calc > t]), 0)),
  median_LFC = sapply(c(5, 4, 3.5, 3, 2.5),
                      function(t) round(median(
                        all_up$mean_abs_LFC[all_up$score_calc > t]), 2))
)
write_csv(dist_table,
          file.path(OUT_DIR, "score_distribution_table.csv"))
cat("\nScore distribution table saved.\n")

# =============================================================================
# SELECT GENES ABOVE THRESHOLD
# =============================================================================
selected <- all_up %>%
  filter(score_calc > SCORE_THRESHOLD)

cat(sprintf("\nSelected genes (score > %.1f): %d\n",
            SCORE_THRESHOLD, nrow(selected)))
cat(sprintf("  baseMean — min: %.0f | median: %.0f | mean: %.0f\n",
            min(selected$mean_baseMean),
            median(selected$mean_baseMean),
            mean(selected$mean_baseMean)))
cat(sprintf("  LFC — min: %.2f | median: %.2f | mean: %.2f\n",
            min(selected$mean_abs_LFC),
            median(selected$mean_abs_LFC),
            mean(selected$mean_abs_LFC)))
cat("  Cluster breakdown:\n")
print(table(selected$cluster))

# =============================================================================
# ANNOTATE WITH GENE NAMES AND NLR/PRR IDENTITY
# =============================================================================
selected_ann <- selected %>%
  left_join(gene_names, by = "gene_id") %>%
  mutate(
    identity = case_when(
      !is.na(nlr_class) & nlr_class != "" ~ "NLR",
      !is.na(prr_class) & prr_class != "" ~ "PRR",
      TRUE ~ "Other"
    )
  )

cat("\nIdentity breakdown:\n")
print(table(selected_ann$identity))

cat("\nNLRs in selected set:\n")
nlrs <- filter(selected_ann, identity == "NLR")
if (nrow(nlrs) > 0) {
  print(nlrs %>% select(gene_id, nlr_class, cluster,
                         mean_baseMean, mean_abs_LFC, n_sig, peak_contrast))
} else {
  cat("  None\n")
}

cat("\nPRRs in selected set:\n")
prrs <- filter(selected_ann, identity == "PRR")
if (nrow(prrs) > 0) {
  print(prrs %>% select(gene_id, prr_class, cluster,
                         mean_baseMean, mean_abs_LFC, n_sig, peak_contrast))
} else {
  cat("  None\n")
}

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

# Full annotated table
write_csv(selected_ann,
          file.path(OUT_DIR, "marker_genes_score3.csv"))

# Gene ID list for promoter extraction
write_lines(selected_ann$gene_id,
            file.path(OUT_DIR, "score3_up_genes.txt"))

# Excel workbook with multiple sheets
sheets <- list(
  "All_selected"     = selected_ann,
  "NLR_PRR_only"     = filter(selected_ann, identity != "Other"),
  "Score_dist_table" = dist_table
)
write_xlsx(sheets,
           file.path(OUT_DIR, "marker_genes_annotated.xlsx"))

cat(sprintf("\nOutputs saved to: %s\n", OUT_DIR))
cat("  marker_genes_score3.csv\n")
cat("  score3_up_genes.txt\n")
cat("  score_distribution_table.csv\n")
cat("  marker_genes_annotated.xlsx\n")
cat("\nDone.\n")
