# Nicotiana benthamiana Immune Receptor RNA-seq Pipeline

![Snakemake](https://img.shields.io/badge/snakemake-≥7.32-brightgreen)
![License](https://img.shields.io/badge/license-MIT-blue)

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
│
├── Reproducible RNA-seq pipeline/     # ── Snakemake pipeline (plug-and-play) ──
│   ├── Snakefile                      #    8-rule workflow: QC→trim→align→count→DESeq2
│   ├── config.yaml                    #    all parameters (fill in genome paths + reference)
│   ├── environment.yml                #    pinned conda environment
│   ├── samples_template.csv           #    copy to samples.csv and fill in
│   ├── scripts/
│   │   ├── deseq2_analysis.R          #    DESeq2 + PCA, volcano, MA, heatmap
│   │   └── qc_plots.R                 #    publication-quality QC bar charts
│   └── example_output/                #    real output from Ma et al. (2025) validation run
│
├── Promoters and HOMER/               # ── Promoter extraction + motif enrichment ──
│   ├── run_homer_all8.sh              #    all 8 HOMER pairwise comparisons (JASPAR2024)
│   ├── extract_promoters.py           #    NLR promoter extraction
│   ├── extract_all_prr_promoters.py   #    PRR promoter extraction
│   ├── extract_subsets.py             #    condition-specific subsets
│   └── extract_53.py
│
├── PRJNA945175/                       # ── FYP viral infection dataset ──
│   ├── align_all.sh                   #    STAR alignment scripts
│   ├── download_and_trim.sh           #    download + Trim Galore
│   ├── featurecounts.sh               #    featureCounts counting
│   ├── nlr_prr_full_landscape.xlsx    #    NLR/PRR classification table
│   ├── nlr_prr_clean_summary.xlsx     #    cleaned NLR/PRR summary
│   ├── NLR_phylogeneti_classification.xlsx  # NLR subclass assignments
│   ├── figures/
│   │   ├── R/                         #    R scripts for all thesis figures
│   │   └── python/                    #    Python scripts for all thesis figures
│   └── tables/                        #    Python scripts for appendix tables
│
├── NLR and PRR extraction/            # ── receptor identification pipeline ──
│   ├── extract_nlr_prr_interpro.py    #    build NLR/PRR landscape from InterProScan + trackers
│   ├── nlr_prr_clean_excel.R          #    generate expressed/DEG summary across 9 contrasts
│   └── run_blastp_prr_vs_arabidopsis.sh #  BLASTp PRRs vs Arabidopsis Swiss-Prot for gene labels
├── deseq2/
│   └── deseq2_clean_rerun.R           #    full DESeq2 analysis (9 contrasts)
├── exploratory/                       #    MEME, GO enrichment, clust, rMATS (not in report)
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

## Snakemake Pipeline (plug-and-play)

The Snakemake pipeline (`Snakefile`) is a standalone, reusable version of this workflow. It takes raw paired-end FASTQs and produces MultiQC reports, sorted BAMs, a gene count matrix, and full DESeq2 output (DEG tables + PCA / volcano / MA / heatmap).

**Validated with:** Ma et al. (2025), PRJNA936199 — 9 samples, 3 conditions, 50 steps, 14,769 DEGs (D36E vs mock), ~95% mapping rate to NbT2T v12.

### Quick start

```bash
cd "Reproducible RNA-seq pipeline/"

# 1. Create conda environment
mamba env create -f environment.yml
conda activate rnaseq-pipeline

# 2. Fill in your samples
cp samples_template.csv samples.csv
# edit samples.csv: one row per replicate (sample_id, condition, fastq_r1, fastq_r2)

# 3. Edit config.yaml — set genome paths and reference_condition

# 4. Dry run to verify
snakemake -n --cores 1

# 5. Run
snakemake --cores 20
```

See [`Reproducible RNA-seq pipeline/example_output/`](Reproducible%20RNA-seq%20pipeline/example_output/) for what the outputs look like when run on a real dataset (Ma et al. 2025, *N. benthamiana* bacterial infection, 9 samples).

### ⚠️ Strandedness — check before running or counts will be wrong

The default is `strandedness: 2` (reverse-stranded, correct for TruSeq/NEBNext). Wrong strandedness produces silently incorrect counts. After the first alignment:

```bash
awk 'NR>4 {s2+=$2; s3+=$3; s4+=$4}
     END {print "unstranded:", s2, "\nforward:", s3, "\nreverse:", s4}' \
    results/alignment/SAMPLE_NAME/ReadsPerGene.out.tab
```

| Output | `strandedness` value |
|--------|---------------------|
| `reverse` >> `forward` | `2` ✓ default |
| `forward` >> `reverse` | `1` |
| `forward` ≈ `reverse` | `0` |

If you need to change it: update `featurecounts: strandedness:` in `config.yaml`, then:

```bash
snakemake --forcerun featurecounts --cores 20
```

### Pipeline outputs

```
results/
├── qc/raw/multiqc_report.html           ← raw read QC
├── qc/trimmed/multiqc_report_trimmed.html
├── alignment/{sample}/
│   ├── Aligned.sortedByCoord.out.bam    ← sorted BAM
│   └── ReadsPerGene.out.tab             ← strandedness check
├── counts/counts_matrix.tsv             ← gene count matrix
└── deseq2/
    ├── DEGs_full.csv / DEGs_significant.csv
    ├── deseq2_summary.txt
    └── plots/  pca_plot.pdf  volcano_plot.pdf  ma_plot.pdf  heatmap_top50.pdf
```

### Troubleshooting

| Symptom | Likely cause | Fix |
|---------|-------------|-----|
| Mapping rate <60% | Wrong genome or STAR index mismatch | Check genome/GTF version match |
| % Assigned very low | Wrong strandedness | Run strandedness check above |
| DEGs = 0 | Wrong `reference_condition` or strandedness | Check config matches samples.csv |
| `MissingInputException` | Wrong FASTQ path | `ls` each path in samples.csv |
| `IncompleteFilesException` | Pipeline interrupted | `--rerun-incomplete` |

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
