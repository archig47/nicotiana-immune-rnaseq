#!/bin/bash
# ============================================================
# STAR Alignment Pipeline — GSE201377
# All 12 samples, per-sample trimmed deletion
# Uses fixed GTF: NbT2T.annot.fixed2.gtf
# Author: fyp5, Kourelis Lab, Imperial College London
# ============================================================

set -euo pipefail

# --- Paths ---
TRIMMED=~/project/GSE201377/trimmed
ALIGNED=~/project/GSE201377/aligned
GENOME=~/project/genome/STAR_index
GTF=~/project/genome/NbT2T.annot.fixed2.gtf
LOG_DIR=~/project/GSE201377/logs/star

mkdir -p "$LOG_DIR"

# --- Disk check helper ---
check_disk() {
    local USED=$(quota -s | awk 'NR==3{gsub(/M/,"",$1); print int($1/1024)}')
    local FREE=$(( 100 - USED ))
    echo "[DISK] ~${USED}GB used, ~${FREE}GB free" | tee -a "$LOG_DIR/pipeline.log"
    if [ "$FREE" -lt 5 ]; then
        echo "[ERROR] Less than 5GB free — stopping pipeline to protect data." | tee -a "$LOG_DIR/pipeline.log"
        exit 1
    fi
}

# --- Alignment + immediate cleanup per sample ---
run_star() {
    local SAMPLE=$1
    local R1=$2
    local R2=$3

    local OUT="$ALIGNED/$SAMPLE"
    mkdir -p "$OUT"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: $SAMPLE" | tee -a "$LOG_DIR/pipeline.log"
    check_disk

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

    # Verify BAM exists and is non-empty before deleting trimmed files
    if [ ! -s "$OUT/Aligned.sortedByCoord.out.bam" ]; then
        echo "[ERROR] BAM missing or empty for $SAMPLE — NOT deleting trimmed files. Stopping." | tee -a "$LOG_DIR/pipeline.log"
        exit 1
    fi

    # Verify count file has gene entries (more than 10 lines = real counts present)
    local COUNT_LINES=$(wc -l < "$OUT/ReadsPerGene.out.tab")
    if [ "$COUNT_LINES" -lt 100 ]; then
        echo "[ERROR] ReadsPerGene.out.tab has only $COUNT_LINES lines for $SAMPLE — GTF mismatch suspected. NOT deleting trimmed files. Stopping." | tee -a "$LOG_DIR/pipeline.log"
        exit 1
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished: $SAMPLE — BAM verified, $COUNT_LINES genes counted." | tee -a "$LOG_DIR/pipeline.log"

    # Delete this sample's trimmed files immediately
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Deleting trimmed files for $SAMPLE..." | tee -a "$LOG_DIR/pipeline.log"
    rm -v "$R1" "$R2" | tee -a "$LOG_DIR/pipeline.log"

    check_disk
    echo "------------------------------------------------------------" | tee -a "$LOG_DIR/pipeline.log"
}

# ============================================================
# BATCH 1 — Pta samples
# ============================================================
echo "==============================" | tee -a "$LOG_DIR/pipeline.log"
echo "BATCH 1: Pta samples" | tee -a "$LOG_DIR/pipeline.log"
echo "==============================" | tee -a "$LOG_DIR/pipeline.log"

run_star "GSM6062243_Pta_rep6" \
    "$TRIMMED/GSM6062243_Pta_rep6_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062243_Pta_rep6_R2_val_2.fq.gz"

run_star "GSM6062244_Pta_rep5" \
    "$TRIMMED/GSM6062244_Pta_rep5_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062244_Pta_rep5_R2_val_2.fq.gz"

run_star "GSM6062245_Pta_rep4" \
    "$TRIMMED/GSM6062245_Pta_rep4_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062245_Pta_rep4_R2_val_2.fq.gz"

run_star "GSM6062249_Pta_rep3" \
    "$TRIMMED/GSM6062249_Pta_rep3_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062249_Pta_rep3_R2_val_2.fq.gz"

run_star "GSM6062250_Pta_rep2" \
    "$TRIMMED/GSM6062250_Pta_rep2_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062250_Pta_rep2_R2_val_2.fq.gz"

run_star "GSM6062251_Pta_rep1" \
    "$TRIMMED/GSM6062251_Pta_rep1_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062251_Pta_rep1_R2_val_2.fq.gz"

# ============================================================
# BATCH 2 — Psy samples
# ============================================================
echo "==============================" | tee -a "$LOG_DIR/pipeline.log"
echo "BATCH 2: Psy samples" | tee -a "$LOG_DIR/pipeline.log"
echo "==============================" | tee -a "$LOG_DIR/pipeline.log"

run_star "GSM6062246_Psy_rep6" \
    "$TRIMMED/GSM6062246_Psy_rep6_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062246_Psy_rep6_R2_val_2.fq.gz"

run_star "GSM6062247_Psy_rep5" \
    "$TRIMMED/GSM6062247_Psy_rep5_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062247_Psy_rep5_R2_val_2.fq.gz"

run_star "GSM6062248_Psy_rep4" \
    "$TRIMMED/GSM6062248_Psy_rep4_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062248_Psy_rep4_R2_val_2.fq.gz"

run_star "GSM6062252_Psy_rep3" \
    "$TRIMMED/GSM6062252_Psy_rep3_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062252_Psy_rep3_R2_val_2.fq.gz"

run_star "GSM6062253_Psy_rep2" \
    "$TRIMMED/GSM6062253_Psy_rep2_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062253_Psy_rep2_R2_val_2.fq.gz"

run_star "GSM6062254_Psy_rep1" \
    "$TRIMMED/GSM6062254_Psy_rep1_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062254_Psy_rep1_R2_val_2.fq.gz"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] ALL 12 SAMPLES COMPLETE." | tee -a "$LOG_DIR/pipeline.log"
check_disk
