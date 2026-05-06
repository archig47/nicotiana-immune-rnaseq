#!/bin/bash
set -e

STAR=~/miniconda3/envs/nlr_pipeline/bin/STAR
TRIMMED=/users/fyp/fyp5/project/GSE201377/trimmed
ALIGNED=/users/fyp/fyp5/project/GSE201377/aligned
INDEX=/users/fyp/fyp5/project/genome/STAR_index_softmasked

mkdir -p $ALIGNED

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

echo "=== STAR alignment started $(date) ==="
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Aligning $SAMPLE..."
    mkdir -p $ALIGNED/$SAMPLE
    $STAR \
        --runThreadN 20 \
        --genomeDir $INDEX \
        --readFilesIn \
            $TRIMMED/${SAMPLE}_R1_val_1.fq.gz \
            $TRIMMED/${SAMPLE}_R2_val_2.fq.gz \
        --readFilesCommand zcat \
        --sjdbOverhang 100 \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMultimapNmax 20 \
        --alignIntronMin 21 \
        --alignIntronMax 1000000 \
        --quantMode GeneCounts \
        --outSJtype Standard \
        --outSAMunmapped None \
        --outReadsUnmapped None \
        --outFileNamePrefix $ALIGNED/$SAMPLE/
    echo "Done $SAMPLE"
done
echo "=== STAR alignment finished $(date) ==="
