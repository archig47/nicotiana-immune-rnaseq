#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

GFFREAD=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/gffread
GFF=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3
GENOME=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.softmasked.genome.fasta
NLRTRACKER=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.NLRtracker_classified.tsv
OUTDIR=/users/fyp/fyp5/project/genome/nlr_promoters

mkdir -p $OUTDIR

# Step 1: Extract NLR gene IDs (strip -mRNA suffix)
echo "Extracting NLR gene IDs..."
tail -n +2 $NLRTRACKER | \
    awk '{print $1}' | \
    sed 's/-mRNA$//' | \
    sort -u > $OUTDIR/nlr_gene_ids.txt
echo "$(wc -l < $OUTDIR/nlr_gene_ids.txt) NLR gene IDs extracted"

# Step 2: Extract NLR gene features from GFF3
echo "Extracting NLR gene records from GFF3..."
grep -F -f $OUTDIR/nlr_gene_ids.txt $GFF | \
    awk '$3 == "gene"' > $OUTDIR/nlr_genes.gff3
echo "$(wc -l < $OUTDIR/nlr_genes.gff3) gene records found"

# Step 3: Extract 2kb upstream sequences using gffread
echo "Extracting 2kb upstream promoter sequences..."
$GFFREAD $OUTDIR/nlr_genes.gff3 \
    -g $GENOME \
    -U \
    -W \
    --upstream-end 2000 \
    -w $OUTDIR/nlr_promoters_2kb.fa

echo "Done. Output: $OUTDIR/nlr_promoters_2kb.fa"
wc -l $OUTDIR/nlr_promoters_2kb.fa
