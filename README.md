# RNA-seq Snakemake Pipeline

## Quick Start

```bash
conda activate rna-seq-pipeline
snakemake results/deseq2/volcano_plot.pdf --cores 20 --use-conda
```

## Reference
- STAR v2.7.11b for alignment
- featureCounts v2.1.1 for counting
- DESeq2 v1.50.2 for differential expression
