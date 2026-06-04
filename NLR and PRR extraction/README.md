# NLR and PRR Extraction

This folder contains scripts used to identify, classify, and annotate NLR and PRR
immune receptor genes in the *Nicotiana benthamiana* NbT2T v12 genome.

---

## Scripts

| Script | Purpose |
|--------|---------|
| `extract_nlr_prr_interpro.py` | Parses InterProScan JSON output alongside NLRtracker and PRRtracker classifications to build the full NLR/PRR receptor landscape Excel file (`nlr_prr_full_landscape.xlsx`). |
| `nlr_prr_clean_excel.R` | Loads the landscape Excel and all 9 DESeq2 contrast results to produce `nlr_prr_clean_summary.xlsx` — tracking which receptors are annotated, expressed, and differentially expressed. |
| `run_blastp_prr_vs_arabidopsis.sh` | Runs BLASTp of 1,252 PRR protein sequences against the *Arabidopsis thaliana* UniProt Swiss-Prot database (16,418 reviewed proteins) to assign closest Arabidopsis homologue labels used in Figure 8. |

---

## Data availability note

NLRtracker and PRRtracker classifications, and NLR phylogenetic subclass assignments,
were performed by the Kourelis laboratory and are not included in this repository.
These tools are described in Kourelis et al. (2021), *eLife*; the v12 annotation
applied here is unpublished. Please contact the Kourelis laboratory directly for
access to these outputs.
