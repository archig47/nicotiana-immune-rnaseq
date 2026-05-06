#!/bin/bash
set -euo pipefail

TRIMMED=~/project/GSE201377/trimmed
ALIGNED=~/project/GSE201377/aligned
GENOME=~/project/genome/STAR_index
GTF=~/project/genome/NbT2T.annot.fixed2.gtf
LOG_DIR=~/project/GSE201377/logs/star

run_star() {
    local SAMPLE=$1
    local R1=$2
    local R2=$3
    local OUT="$ALIGNED/$SAMPLE"
    mkdir -p "$OUT"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: $SAMPLE" | tee -a "$LOG_DIR/pipeline.log"

    STAR \
        --runThreadN 20 \
        --genomeDir "$GENOME" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang 100 \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFilterMultimapNmax 20 \
        --alignIntronMin 21 \
        --alignIntronMax 1000000 \
        --quantMode GeneCounts \
        --outSJtype Standard \
        --outSAMunmapped None \
        --outReadsUnmapped None \
        --outBAMsortingThreadN 0 \
        --outFileNamePrefix "$OUT/" \
        2>&1 | tee "$LOG_DIR/${SAMPLE}.log"

    if [ ! -s "$OUT/Aligned.sortedByCoord.out.bam" ]; then
        echo "[ERROR] BAM missing or empty for $SAMPLE — stopping." | tee -a "$LOG_DIR/pipeline.log"
        exit 1
    fi

    local COUNT_LINES=$(wc -l < "$OUT/ReadsPerGene.out.tab")
    if [ "$COUNT_LINES" -lt 100 ]; then
        echo "[ERROR] Only $COUNT_LINES lines in count file for $SAMPLE — stopping." | tee -a "$LOG_DIR/pipeline.log"
        exit 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished: $SAMPLE — BAM verified, $COUNT_LINES genes counted." | tee -a "$LOG_DIR/pipeline.log"
    rm -v "$R1" "$R2" | tee -a "$LOG_DIR/pipeline.log"
}

run_star "GSM6062253_Psy_rep2" \
    "$TRIMMED/GSM6062253_Psy_rep2_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062253_Psy_rep2_R2_val_2.fq.gz"

run_star "GSM6062254_Psy_rep1" \
    "$TRIMMED/GSM6062254_Psy_rep1_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062254_Psy_rep1_R2_val_2.fq.gz"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ALL 3 REMAINING SAMPLES COMPLETE." | tee -a "$LOG_DIR/pipeline.log"
