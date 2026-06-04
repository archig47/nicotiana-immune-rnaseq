# Exploratory Analyses

These analyses were conducted during the course of the FYP project to investigate
the *N. benthamiana* immune receptor transcriptome more broadly. They were not
included in the final report but were instrumental in exploring the biology,
testing analytical approaches, and understanding the genome and dataset.

They are preserved here for completeness and as a resource for future work.

---

## Contents

### `meme/` — MEME Suite motif analysis (AME + STREME)
AME (Analysis of Motif Enrichment) and STREME (de novo motif discovery) were run
on PRR promoter subsets using the MEME Suite v5.5.8 prior to the decision to use
HOMER with the JASPAR2024 library. AME results are provided for condition-specific
PRR subsets (HCRV up/down, TSWV up/down, TH up/down, pan-viral up). STREME de
novo results are also included. HOMER was ultimately selected for the final
analysis due to its direct support for pairwise foreground/background comparisons
(see Appendix 6.2.1 of the thesis for full rationale).

### `clust/` — Marker gene selection
`select_marker_genes.R` implements a composite score-based approach for ranking
DEGs by response magnitude and breadth, used in exploratory prioritisation of
candidate NLR and PRR genes for functional follow-up.

### `go_enrichment/` — Gene Ontology enrichment
GO enrichment was performed on condition-specific DEG subsets (HCRV upregulated,
TSWV downregulated, TSWV-specific downregulated) using both agriGO and gProfiler2.
Results and dot plot visualisations are included. These analyses confirmed the
expected immune and metabolic pathway enrichments but were not incorporated into
the final report due to space constraints.

### `rmats/` — Alternative splicing analysis (rMATS)
rMATS v4.1 was used to detect differential alternative splicing in NLR and PRR
genes across all nine viral infection conditions. Scripts for parsing results,
and generating summary heatmaps and dot plots are included. This analysis was
completed but removed from the final report scope; NbT2T18g02284 (significant in
all 9 conditions) was identified as a strong candidate for future functional
validation.

---

*None of these analyses appear in the final FYP report. They are included here
as a record of exploratory work and may serve as a starting point for future
investigations.*
