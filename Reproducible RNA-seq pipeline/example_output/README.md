# Example Output — Ma et al. (2025), PRJNA936199

This folder shows what the pipeline produces when run on a real dataset.

**Dataset:** *Nicotiana benthamiana* bacterial infection RNA-seq  
**Reference:** Ma et al. (2025), BioProject PRJNA936199  
**Conditions:** mock (3 reps) | D36E — flagellin+ *P. syringae* (3 reps) | dfliC — flagellin− mutant (3 reps)  
**Genome:** NbT2T v12

---

## Files in this folder

| File | What it is |
|------|-----------|
| `deseq2_summary.txt` | Plain-text summary: genes tested, DEG counts per contrast, QC table |
| `plots/pca_plot.png` | PCA of VST-normalised counts — conditions separate cleanly on PC1 (91% variance) |
| `plots/volcano_D36E_vs_mock.png` | Volcano plot — 14,769 DEGs in D36E vs mock |
| `plots/heatmap_top50.png` | Top 50 DEGs by adjusted p-value, clustered |
| `plots/qc_read_counts.png` | Raw read depth per sample (all >20M) |
| `plots/qc_mapping_rate.png` | Unique mapping rate per sample (all ~95%) |
| `plots/qc_assigned_reads.png` | % reads assigned to genes (all ~83%) |
| `plots/supplementary_figure_pipeline_validation.png` | 6-panel supplementary figure (A–F) |

---

## Key results

- **36,554 genes** tested after low-count filter
- **D36E vs mock:** 7,728 up + 7,041 down = **14,769 DEGs**  
  (D36E has functional flagellin → strong PTI response)
- **dfliC vs mock:** 2,291 up + 655 down = **2,946 DEGs**  
  (flagellin− mutant → attenuated response, confirms flagellin-dependent signalling)
- All 9 samples pass QC: >20M reads, ~95% uniquely mapped, ~83% assigned to genes

---

## How to reproduce

```bash
cd pipeline/
cp samples_template.csv samples.csv
# fill in samples.csv with your FASTQ paths
# edit config.yaml: set genome paths, reference_condition: "mock"
snakemake --cores 20
```

See the main [pipeline README](../README.md) and the top-level [repository README](../../README.md) for full instructions.
