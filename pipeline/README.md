# RNA-seq Pipeline

Plug-and-play Snakemake pipeline: paired-end FASTQs → MultiQC → STAR alignment → featureCounts → DESeq2.

See [`example_output/`](example_output/) for real outputs from the validation run (Ma et al. 2025, PRJNA936199).

---

## Quick start

```bash
# 1. Install environment (once)
mamba env create -f environment.yml
conda activate rnaseq-pipeline

# 2. Fill in samples
cp samples_template.csv samples.csv
# edit samples.csv — one row per replicate

# 3. Set genome paths + reference condition in config.yaml

# 4. Dry run
snakemake -n --cores 1

# 5. Run
snakemake --cores 20
```

## What you get

```
results/
├── qc/raw/multiqc_report.html
├── qc/trimmed/multiqc_report_trimmed.html
├── alignment/{sample}/
│   ├── Aligned.sortedByCoord.out.bam
│   └── ReadsPerGene.out.tab        ← use for strandedness check
├── counts/counts_matrix.tsv
└── deseq2/
    ├── DEGs_full.csv
    ├── DEGs_significant.csv
    ├── deseq2_summary.txt
    └── plots/
        ├── pca_plot.pdf
        ├── volcano_plot.pdf
        ├── ma_plot.pdf
        └── heatmap_top50.pdf
```

## ⚠️ Strandedness

The default `strandedness: 2` (reverse-stranded) is correct for TruSeq/NEBNext kits. Wrong strandedness silently produces incorrect counts. After the first alignment, verify:

```bash
awk 'NR>4 {s2+=$2; s3+=$3; s4+=$4}
     END {print "unstranded:", s2, "\nforward:", s3, "\nreverse:", s4}' \
    results/alignment/SAMPLE/ReadsPerGene.out.tab
```

| Output | `strandedness` |
|--------|---------------|
| reverse >> forward | `2` ✓ default |
| forward >> reverse | `1` |
| forward ≈ reverse | `0` |

Change in `config.yaml` then re-run: `snakemake --forcerun featurecounts --cores 20`

## Troubleshooting

| Symptom | Fix |
|---------|-----|
| Mapping rate <60% | Check genome/GTF version match |
| % Assigned very low | Check strandedness |
| DEGs = 0 | Check `reference_condition` in config matches samples.csv |
| `MissingInputException` | `ls` each FASTQ path in samples.csv |
| `IncompleteFilesException` | `snakemake --rerun-incomplete --cores 20` |

Full documentation in the [repository README](../README.md).
