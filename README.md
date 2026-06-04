# Nicotiana benthamiana Immune Receptor RNA-seq Pipeline

**Author:** Archita Gupta  
**Lab:** Kourelis Lab, Imperial College London  
**Project:** FYP 2025–26

---

## Overview

Reproducible RNA-seq pipeline characterising the transcriptional and cis-regulatory landscape of NLR and PRR immune receptors in *Nicotiana benthamiana* during viral infection. The pipeline was developed and validated against the Kourelis lab GFF3 annotation (v12) of the NbT2T telomere-to-telomere genome assembly.

---

## Datasets

| Accession | Description | Conditions |
|-----------|-------------|------------|
| PRJNA945175 | Viral infection time course | HCRV, TSWV, TH co-infection vs mock (CK) at 1, 7, 14 dpi |
| PRJNA936199 | Pipeline validation (Ma et al., 2025) | *P. syringae* DC3000 D36E vs mock |

---

## Pipeline Steps

1. **Download & QC** — `fasterq-dump` (SRA Toolkit v3.2.1), Falco v1.2.5, MultiQC v1.33
2. **Adapter trimming** — Trim Galore v0.6.11 (Q20, 30 bp minimum, stringency 3)
3. **Alignment** — STAR v2.7.11b against NbT2T genome (v12 annotation)
4. **Read counting** — featureCounts v2.1.1 (reverse-stranded, fractional multimapper counting)
5. **Differential expression** — DESeq2 v1.50.2 in R v4.5.3 (9 independent virus-vs-mock contrasts; padj < 0.05, |log₂FC| > 1)
6. **NLR/PRR classification** — NLRtracker/PRRtracker (Kourelis et al., 2021); NLR subclasses by maximum-likelihood phylogeny (NB-ARC domains, ClustalO + FastTree in Geneious)
7. **Promoter extraction** — BEDTools v2.31.1 (2 kb upstream TSS, strand-aware)
8. **Motif enrichment** — HOMER v5.1 with JASPAR2024 CORE plants library (805 motifs); 8 pairwise foreground/background comparisons

---

## Repository Structure

```
nicotiana-immune-rnaseq/
├── PRJNA945175/              # Shell scripts for download, alignment, counting
├── genome/                   # Genome-related resources
├── deseq2/
│   └── deseq2_clean_rerun.R  # Full DESeq2 analysis script
├── figures/
│   ├── R/                    # R scripts for figure generation
│   │   ├── pca_final.R
│   │   ├── prr_heatmap_split_v3.R
│   │   ├── prr_upset_v5.R
│   │   └── nlr_prr_figures_v2.R
│   └── python/               # Python scripts for figure generation
│       ├── prr_dotplot_v10.py
│       ├── nlr_scatter_v5.py
│       ├── prr_scatter_v1.py
│       ├── homer_dotplot_v6.py
│       └── qc_figures.py
└── README.md
```

---

## Key Results

- **666 PRR DEGs** and **98 NLR DEGs** identified across 9 virus-vs-mock contrasts
- PRR DEG rate significantly enriched vs expressed genome background (OR = 1.46, p = 1.57 × 10⁻¹⁴)
- WRKY motifs dominant across both receptor families; MYB motifs NLR-specific; DREB1/CBF and DOF4.2 PRR-associated
- TH co-infection produced attenuated responses relative to either single infection
- CC-clade2 NLRs showed strongest upregulation; LRR-RLKs showed TSWV-specific temporal inversion

---

## Reference Genome

- **Assembly:** NbT2T (Chen et al., 2024, *Nature Plants*)
- **Annotation:** Kourelis lab GFF3 v12 (99% BUSCO completeness, unpublished)
- **NLR/PRR inventory:** 284 NLRs, 1,252 PRRs identified by NLRtracker/PRRtracker

---

## InterProScan

Domain annotations were generated using InterProScan v5.77-108.0 against the full NbT2T v12 proteome. To regenerate:

```bash
interproscan.sh -i proteins_v12.faa -f tsv -dp -goterms --cpu 20
```

---

## Requirements

- STAR v2.7.11b
- featureCounts (Subread v2.1.1)
- Trim Galore v0.6.11
- Falco v1.2.5
- DESeq2 v1.50.2 (R v4.5.3)
- HOMER v5.1
- BEDTools v2.31.1
- Python 3 (matplotlib, pandas, numpy, scipy, logomaker, adjustText)
- conda env: `nlr_pipeline`

---

## License

MIT License
