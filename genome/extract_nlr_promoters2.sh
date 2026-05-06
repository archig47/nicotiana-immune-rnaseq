#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

GENOME=/users/fyp/fyp5/project/genome/jiorgos/NbT2T.softmasked.genome.fasta
OUTDIR=/users/fyp/fyp5/project/genome/nlr_promoters
IDFILE=$OUTDIR/nlr_gene_ids.txt

echo "Creating BED file of 2kb upstream regions..."

# Parse GFF3 gene records for NLR genes, extract 2kb upstream (strand-aware)
grep -F -f $IDFILE /users/fyp/fyp5/project/genome/jiorgos/NbT2T.final_v12.gff3 | \
awk '$3 == "gene"' | \
awk 'BEGIN{OFS="\t"} {
    # Get gene ID from attributes
    match($9, /ID=([^;]+)/, arr)
    id = arr[1]
    chr = $1
    start = $4 - 1  # convert to 0-based
    end = $5
    strand = $7

    if (strand == "+") {
        # upstream = 2kb before start
        up_start = start - 2000
        if (up_start < 0) up_start = 0
        up_end = start
    } else {
        # upstream = 2kb after end (reverse strand)
        up_start = end
        up_end = end + 2000
    }
    print chr, up_start, up_end, id, ".", strand
}' > $OUTDIR/nlr_promoters_2kb.bed

echo "$(wc -l < $OUTDIR/nlr_promoters_2kb.bed) promoter regions in BED file"

# Extract sequences using samtools faidx regions
echo "Indexing genome..."
samtools faidx $GENOME

echo "Extracting sequences..."
> $OUTDIR/nlr_promoters_2kb.fa

while IFS=$'\t' read -r chr start end name score strand; do
    # samtools faidx uses 1-based coordinates
    region="${chr}:$((start+1))-${end}"
    seq=$(samtools faidx $GENOME "$region" 2>/dev/null | tail -n +2 | tr -d '\n')
    
    if [ -z "$seq" ]; then
        echo "WARNING: No sequence for $name ($region)" >&2
        continue
    fi
    
    # Reverse complement if minus strand
    if [ "$strand" == "-" ]; then
        seq=$(echo "$seq" | tr 'ACGTacgt' 'TGCAtgca' | rev)
    fi
    
    echo ">$name" >> $OUTDIR/nlr_promoters_2kb.fa
    echo "$seq" >> $OUTDIR/nlr_promoters_2kb.fa
done < $OUTDIR/nlr_promoters_2kb.bed

echo "Done."
grep -c "^>" $OUTDIR/nlr_promoters_2kb.fa
