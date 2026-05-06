# Nicotiana benthamiana Immune Receptor RNA-seq Pipeline

**Author:** Archita Gupta  
**Lab:** Kourelis Lab, Imperial College London  
**Project:** FYP 2025-26

## Overview

Reproducible RNA-seq pipeline characterising the transcriptional and cis-regulatory landscape of NLR and PRR immune receptors in *Nicotiana benthamiana* across bacterial and viral infection contexts.

## Datasets

| Dataset | Description | Conditions |
|---------|-------------|------------|
| GSE201377 | Bacterial infection | Psy B728a vs Pta DC3000, no mock |
| PRJNA945175 | Viral infection | HCRV, TSWV, TH co-infection vs CK mock |

## Pipeline Steps

1. Download & QC (TrimGalore, Falco, MultiQC)
2. Alignment (STAR, NbT2T v12 genome)
3. Feature counting (featureCounts)
4. Differential expression (DESeq2, padj < 0.05, |LFC| > 1)
5. NLR/PRR filtering (Kourelis lab NLRtracker/PRRtracker)
6. Promoter extraction (2kb upstream TSS)
7. Motif enrichment (MEME Suite AME + STREME, Arabidopsis DAPv1)

## Key Findings

- HCRV-upregulated PRRs enriched for WRKY motifs (WRKY59, WRKY26, WRKY28) — immune activation signal
- TSWV-specifically suppressed PRRs enriched for bZIP/ABA-responsive motifs (ABI5, ABF2, AREB3) — consistent with active immune evasion
- Null result across full PRR/NLR backgrounds is meaningful: transcriptional selectivity is detectable only within coherent directional subsets

## Genome

NbT2T v12 (Kourelis lab annotation). Softmasked genome used for alignment.  
Motif database: Arabidopsis DAP-seq v1 (`ArabidopsisDAPv1.meme`)

## Requirements

- STAR 2.7+
- featureCounts (subread)
- TrimGalore
- DESeq2 (R)
- MEME Suite 5.5.8
- Python 3, conda env: `nlr_pipeline`
