#!/bin/bash
module load meme/5.5.8

SUBSETS=/users/fyp/fyp5/project/PRJNA945175/promoters/subsets
BG=/users/fyp/fyp5/project/PRJNA945175/promoters/all_PRR_promoters_2kb.fa
MOTIFS=/users/fyp/fyp5/project/genome/motif_databases/ARABD/ArabidopsisDAPv1.meme
OUTBASE=/users/fyp/fyp5/project/PRJNA945175/meme/subsets

mkdir -p $OUTBASE

for fa in $SUBSETS/*.fa; do
    name=$(basename $fa _promoters_2kb.fa)
    outdir=$OUTBASE/${name}_AME
    mkdir -p $outdir
    echo "Running AME on $name..."
    ame \
        --oc $outdir \
        --scoring avg \
        --method fisher \
        --hit-lo-fraction 0.25 \
        --evalue-report-threshold 10.0 \
        --control $BG \
        --kmer 2 \
        $fa \
        $MOTIFS \
    > $outdir/ame.log 2>&1
    hits=$(grep -v "^#" $outdir/ame.tsv | grep -v "^rank" | wc -l)
    echo "  $name: $hits significant motifs"
done

echo "All done."

