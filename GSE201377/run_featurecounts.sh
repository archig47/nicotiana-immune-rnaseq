#!/bin/bash
set -e

FEATURECOUNTS=~/miniconda3/envs/nlr_pipeline/bin/featureCounts
ALIGNED=/users/fyp/fyp5/project/GSE201377/aligned
OUTDIR=/users/fyp/fyp5/project/GSE201377/featurecounts_softmasked
GFF_V12=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3
GFF_HELIXER=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.helixer.gff3

mkdir -p $OUTDIR

BAMS=$(ls $ALIGNED/*/Aligned.sortedByCoord.out.bam | tr '\n' ' ')

echo "=== featureCounts v12 started $(date) ==="
$FEATURECOUNTS \
    -T 12 -p -s 2 \
    -t exon -g Parent \
    -a $GFF_V12 \
    -F GFF \
    -o $OUTDIR/counts_matrix_v12.txt \
    $BAMS
echo "=== featureCounts v12 done $(date) ==="

echo "=== featureCounts v12 multimapper started $(date) ==="
$FEATURECOUNTS \
    -T 12 -p -s 2 \
    -t exon -g Parent \
    -M --fraction \
    -a $GFF_V12 \
    -F GFF \
    -o $OUTDIR/counts_matrix_v12_multimapper.txt \
    $BAMS
echo "=== featureCounts v12 multimapper done $(date) ==="

echo "=== featureCounts Helixer started $(date) ==="
$FEATURECOUNTS \
    -T 12 -p -s 2 \
    -t exon -g Parent \
    -a $GFF_HELIXER \
    -F GFF \
    -o $OUTDIR/counts_matrix_helixer.txt \
    $BAMS
echo "=== featureCounts Helixer done $(date) ==="

echo "=== All featureCounts done $(date) ==="
