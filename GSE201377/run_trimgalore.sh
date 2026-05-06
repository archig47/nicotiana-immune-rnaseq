#!/bin/bash
set -e

TRIMGALORE=~/miniconda3/envs/nlr_pipeline/bin/trim_galore
CUTADAPT=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/cutadapt
MERGED=/users/fyp/fyp5/project/GSE201377/merged
TRIMMED=/users/fyp/fyp5/project/GSE201377/trimmed

mkdir -p $TRIMMED

SAMPLES=(
GSM6062243_Pta_rep6
GSM6062244_Pta_rep5
GSM6062245_Pta_rep4
GSM6062246_Psy_rep6
GSM6062247_Psy_rep5
GSM6062248_Psy_rep4
GSM6062249_Pta_rep3
GSM6062250_Pta_rep2
GSM6062251_Pta_rep1
GSM6062252_Psy_rep3
GSM6062253_Psy_rep2
GSM6062254_Psy_rep1
)

echo "=== Trim Galore started $(date) ==="
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Trimming $SAMPLE..."
    $TRIMGALORE --paired --quality 20 --length 36 --fastqc --cores 4 \
        --path_to_cutadapt $CUTADAPT \
        --output_dir $TRIMMED \
        $MERGED/${SAMPLE}_R1.fastq.gz \
        $MERGED/${SAMPLE}_R2.fastq.gz
    echo "Done $SAMPLE"
done
echo "=== Trim Galore finished $(date) ==="
