#!/bin/bash
set -e

PREFETCH=~/miniconda3/envs/nlr_pipeline/bin/prefetch
FASTERQ=~/miniconda3/envs/nlr_pipeline/bin/fasterq-dump
TRIMGALORE=~/miniconda3/envs/nlr_pipeline/bin/trim_galore
CUTADAPT=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/cutadapt
TMPRAW=/users/fyp/fyp5/project/PRJNA945175/raw_tmp
TRIMMED=/users/fyp/fyp5/project/PRJNA945175/trimmed

mkdir -p $TMPRAW

SAMPLES=(
SRR23875052
SRR23875053
SRR23875054
SRR23875055
SRR23875056
SRR23875057
SRR23875058
SRR23875060
SRR23875061
SRR23875062
SRR23875063
SRR23875064
SRR23875065
SRR23875066
SRR23875067
SRR23875068
SRR23875069
SRR23875070
SRR23875072
SRR23875073
SRR23875074
SRR23875075
SRR23875077
SRR23875079
SRR23875080
SRR23875081
SRR23875082
SRR23875083
)

echo "=== Download and trim started $(date) ==="
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Downloading $SAMPLE..."
    $PREFETCH --max-size 50G $SAMPLE --output-directory $TMPRAW
    $FASTERQ --split-files --threads 8 --outdir $TMPRAW $TMPRAW/$SAMPLE/$SAMPLE.sra
    gzip $TMPRAW/${SAMPLE}_1.fastq
    gzip $TMPRAW/${SAMPLE}_2.fastq
    rm -rf $TMPRAW/$SAMPLE

    echo "Trimming $SAMPLE..."
    $TRIMGALORE --paired --quality 20 --length 36 --fastqc --cores 4 \
        --path_to_cutadapt $CUTADAPT \
        --output_dir $TRIMMED \
        $TMPRAW/${SAMPLE}_1.fastq.gz \
        $TMPRAW/${SAMPLE}_2.fastq.gz

    rm $TMPRAW/${SAMPLE}_1.fastq.gz $TMPRAW/${SAMPLE}_2.fastq.gz
    echo "Done $SAMPLE"
done

rmdir $TMPRAW
echo "=== All done $(date) ==="
