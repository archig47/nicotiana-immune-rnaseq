# PRJNA945175 — Viral Infection Time Course

This directory contains all scripts, input data, and analysis code for the primary
FYP dataset: a publicly available *Nicotiana benthamiana* RNA-seq time course of
viral infection (Gui et al., 2023; NCBI BioProject PRJNA945175).

**Conditions:** HCRV, TSWV, and TH co-infection vs mock (CK) at 1, 7, and 14 dpi
**Samples:** 30 libraries (2–3 biological replicates per condition–timepoint combination)

---

## NLR/PRR Classification Scripts

These scripts are located in `../NLR and PRR extraction/` at the root of the repository.

| Script | Purpose |
|--------|---------|
| `extract_nlr_prr_interpro.py` | Parses InterProScan JSON output alongside NLRtracker and PRRtracker classifications to build the full NLR/PRR receptor landscape Excel file (`nlr_prr_full_landscape.xlsx`). |
| `nlr_prr_clean_excel.R` | Loads the landscape Excel and all 9 DESeq2 contrast results to produce `nlr_prr_clean_summary.xlsx`. |
| `run_blastp_prr_vs_arabidopsis.sh` | Runs BLASTp of 1,252 PRR protein sequences against the *Arabidopsis thaliana* UniProt Swiss-Prot database to assign gene labels used in Figure 8. |

---

## Pipeline Scripts

| Script | Purpose |
|--------|---------|
| `download_and_trim.sh` | Download raw reads via fasterq-dump and run Trim Galore |
| `align_all.sh` | Align all 30 samples to NbT2T genome using STAR |
| `align_test.sh` | Test alignment on a single sample |
| `align_topup.sh` | Re-align samples that failed or were added later |
| `featurecounts.sh` | Count reads per gene using featureCounts (fractional multimapper mode) |
| `recover_081.sh` | Recovery script for SRR23875081 (processed separately) |

The full DESeq2 differential expression analysis is in `../deseq2/deseq2_clean_rerun.R`,
which produces 9 independent virus-vs-mock contrasts (padj < 0.05, \|log₂FC\| > 1).

---

## Input Data

| File | Contents |
|------|---------|
| `nlr_prr_full_landscape.xlsx` | Full NLR and PRR receptor classification (ectodomain, backbone, matched rules) for all 1,536 annotated receptors in the NbT2T v12 proteome |
| `nlr_prr_clean_summary.xlsx` | Cleaned summary of expressed and differentially expressed NLRs and PRRs |
| `NLR_phylogeneti_classification.xlsx` | NLR subclass assignments derived from maximum-likelihood phylogenetic analysis of NB-ARC domain sequences (CC-clade1/2, TIR, NRC-sensor-01/02, NRC-helper, TNP, Unclassified) |

---

## Figure Scripts

Scripts that generate all figures in the FYP report.

### `figures/R/`

| Script | Figure |
|--------|--------|
| `pca_final.R` | PCA of rlog-normalised counts across all timepoints (Figure 3) |
| `prr_heatmap_split_v3.R` | PRR class-level heatmap split by direction (Figure 7) |
| `prr_upset_v5.R` | PRR DEG overlap UpSet plot (Figure 9) |
| `nlr_upset_custom.R` | NLR DEG overlap UpSet plot (Figure 6) |
| `volcano_all9_capped60.R` | Volcano plots for all 9 contrasts, y-axis capped at −log₁₀p = 60 (Figure 4) |

### `figures/python/`

| Script | Figure |
|--------|--------|
| `nlr_heatmap_v11.py` | NLR DEG expression heatmap by subclass (Figure 5) |
| `nlr_scatter_v5.py` | NLR DEG abundance vs fold change scatter (Appendix Figure D) |
| `prr_dotplot_v10.py` | Top PRR DEGs dotplot by ectodomain class (Figure 8) |
| `prr_scatter_v1.py` | PRR DEG abundance vs fold change scatter (Appendix Figure E) |
| `homer_occurrence_v1.py` | HOMER motif prevalence heatmap (Figure 10, panel I) |
| `homer_position_v2.py` | HOMER motif positional KDE distribution (Figure 10, panels IIA–B) |
| `qc_figures.py` | QC bar charts: read depth, alignment rates, featureCounts assignment (Appendix Figures A–C) |

---

## Appendix Table Scripts

### `tables/`

| Script | Table |
|--------|-------|
| `make_package_tables.py` | R and Python package tables (Tables A.1 and A.2) |
| `make_tool_comparison_table.py` | Tool selection rationale table (Table A.3) |
| `make_homer_design_table.py` | HOMER comparison design table (Table A.4) |
