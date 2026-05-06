#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

GFF_V12=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3
GFF_HELIXER=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.helixer.gff3
OUTDIR=/users/fyp/fyp5/project/PRJNA945175/featurecounts
LOGDIR=/users/fyp/fyp5/project/PRJNA945175/logs
FEATURECOUNTS=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/featureCounts

mkdir -p $OUTDIR

# Build BAM list — all 31 current good BAMs (081 will be added later)
BAMS=$(ls /users/fyp/fyp5/project/PRJNA945175/aligned/*/Aligned.sortedByCoord.out.bam | sort)

echo "=== featureCounts started $(date) ===" | tee -a $LOGDIR/featurecounts.log
echo "BAMs included:" | tee -a $LOGDIR/featurecounts.log
echo "$BAMS" | tee -a $LOGDIR/featurecounts.log

# Run 1 — v12 standard
echo "[$(date)] Starting Run 1: v12 standard..." | tee -a $LOGDIR/featurecounts.log
$FEATURECOUNTS \
    -T 20 \
    -p \
    -s 2 \
    -t exon \
    -g Parent \
    -F GFF \
    -a $GFF_V12 \
    -o $OUTDIR/counts_v12_standard.txt \
    $BAMS >> $LOGDIR/featurecounts.log 2>&1

if [ $? -eq 0 ]; then
    echo "[$(date)] Run 1 complete" | tee -a $LOGDIR/featurecounts.log
else
    echo "[$(date)] FAILED: Run 1" | tee -a $LOGDIR/featurecounts.log
fi

# Run 2 — v12 multimapper (PRIMARY for DESeq2)
echo "[$(date)] Starting Run 2: v12 multimapper..." | tee -a $LOGDIR/featurecounts.log
$FEATURECOUNTS \
    -T 20 \
    -p \
    -s 2 \
    -t exon \
    -g Parent \
    -F GFF \
    -M --fraction \
    -a $GFF_V12 \
    -o $OUTDIR/counts_v12_multimapper.txt \
    $BAMS >> $LOGDIR/featurecounts.log 2>&1

if [ $? -eq 0 ]; then
    echo "[$(date)] Run 2 complete" | tee -a $LOGDIR/featurecounts.log
else
    echo "[$(date)] FAILED: Run 2" | tee -a $LOGDIR/featurecounts.log
fi

# Run 3 — Helixer standard
echo "[$(date)] Starting Run 3: Helixer standard..." | tee -a $LOGDIR/featurecounts.log
$FEATURECOUNTS \
    -T 20 \
    -p \
    -s 2 \
    -t exon \
    -g Parent \
    -F GFF \
    -a $GFF_HELIXER \
    -o $OUTDIR/counts_helixer_standard.txt \
    $BAMS >> $LOGDIR/featurecounts.log 2>&1

if [ $? -eq 0 ]; then
    echo "[$(date)] Run 3 complete" | tee -a $LOGDIR/featurecounts.log
else
    echo "[$(date)] FAILED: Run 3" | tee -a $LOGDIR/featurecounts.log
fi

echo "=== featureCounts complete $(date) ===" | tee -a $LOGDIR/featurecounts.log
