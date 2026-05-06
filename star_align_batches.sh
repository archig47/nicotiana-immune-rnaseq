#!/bin/bash
# ============================================================
# STAR Alignment Pipeline — GSE201377
# Batches 1 (Pta) and 2 (Psy)
# Author: fyp5, Kourelis Lab, Imperial College London
# ============================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# --- Paths ---
TRIMMED=~/project/GSE201377/trimmed
ALIGNED=~/project/GSE201377/aligned
GENOME=~/project/genome/STAR_index
GTF=~/project/genome/NbT2T.annot.gtf
LOG_DIR=~/project/GSE201377/logs/star

mkdir -p "$LOG_DIR"

# --- STAR parameters ---
THREADS=20
OVERHANG=100

run_star() {
    local SAMPLE=$1
    local R1=$2
    local R2=$3

    local OUT="$ALIGNED/$SAMPLE"
    mkdir -p "$OUT"

    echo "[$( date '+%Y-%m-%d %H:%M:%S' )] Starting: $SAMPLE" | tee -a "$LOG_DIR/pipeline.log"

    STAR \
        --runThreadN $THREADS \
        --genomeDir "$GENOME" \
        --sjdbGTFfile "$GTF" \
        --sjdbOverhang $OVERHANG \
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

    echo "[$( date '+%Y-%m-%d %H:%M:%S' )] Finished: $SAMPLE" | tee -a "$LOG_DIR/pipeline.log"
}

# ============================================================
# BATCH 1 — Pta samples
# ============================================================
echo "==============================" | tee -a "$LOG_DIR/pipeline.log"
echo "BATCH 1: Pta samples starting" | tee -a "$LOG_DIR/pipeline.log"
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

# --- Delete Batch 1 trimmed files to free disk space ---
echo "[$( date '+%Y-%m-%d %H:%M:%S' )] Batch 1 complete. Deleting Pta trimmed files..." | tee -a "$LOG_DIR/pipeline.log"

rm -v \
    "$TRIMMED/GSM6062243_Pta_rep6_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062243_Pta_rep6_R2_val_2.fq.gz" \
    "$TRIMMED/GSM6062244_Pta_rep5_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062244_Pta_rep5_R2_val_2.fq.gz" \
    "$TRIMMED/GSM6062245_Pta_rep4_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062245_Pta_rep4_R2_val_2.fq.gz" \
    "$TRIMMED/GSM6062249_Pta_rep3_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062249_Pta_rep3_R2_val_2.fq.gz" \
    "$TRIMMED/GSM6062250_Pta_rep2_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062250_Pta_rep2_R2_val_2.fq.gz" \
    "$TRIMMED/GSM6062251_Pta_rep1_R1_val_1.fq.gz" \
    "$TRIMMED/GSM6062251_Pta_rep1_R2_val_2.fq.gz"

echo "[$( date '+%Y-%m-%d %H:%M:%S' )] Pta trimmed files deleted." | tee -a "$LOG_DIR/pipeline.log"

# ============================================================
# BATCH 2 — Psy samples
# ============================================================
echo "==============================" | tee -a "$LOG_DIR/pipeline.log"
echo "BATCH 2: Psy samples starting" | tee -a "$LOG_DIR/pipeline.log"
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

echo "[$( date '+%Y-%m-%d %H:%M:%S' )] ALL BATCHES COMPLETE." | tee -a "$LOG_DIR/pipeline.log"
echo "BAM files in: $ALIGNED" | tee -a "$LOG_DIR/pipeline.log"
