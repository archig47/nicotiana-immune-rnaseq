#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

TRIMDIR=/users/fyp/fyp5/project/PRJNA945175/trimmed
OUTDIR=/users/fyp/fyp5/project/PRJNA945175/aligned
LOGDIR=/users/fyp/fyp5/project/PRJNA945175/logs
INDEX=/users/fyp/fyp5/project/genome/STAR_index_softmasked
STAR=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/STAR-avx2

SAMPLES=(
SRR23875051 SRR23875052 SRR23875053 SRR23875054 SRR23875055 SRR23875056
SRR23875057 SRR23875058 SRR23875059 SRR23875060 SRR23875061 SRR23875062
SRR23875063 SRR23875064 SRR23875065 SRR23875066 SRR23875067 SRR23875068
SRR23875069 SRR23875070 SRR23875071 SRR23875072 SRR23875073 SRR23875074
SRR23875075 SRR23875076 SRR23875077 SRR23875078 SRR23875079 SRR23875080
SRR23875081 SRR23875082
)

echo "=== PRJNA945175 STAR alignment started $(date) ===" | tee -a $LOGDIR/align_all.log

for SAMPLE in "${SAMPLES[@]}"; do
    # Skip SRR23875051 — already aligned as test
    if [ "$SAMPLE" == "SRR23875051" ]; then
        echo "[$(date)] Skipping $SAMPLE — already aligned" | tee -a $LOGDIR/align_all.log
        continue
    fi

    R1=$TRIMDIR/${SAMPLE}_1_val_1.fq.gz
    R2=$TRIMDIR/${SAMPLE}_2_val_2.fq.gz

    # Check trimmed files exist
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "[$(date)] ERROR: Missing trimmed files for $SAMPLE — skipping" | tee -a $LOGDIR/align_all.log
        continue
    fi

    mkdir -p $OUTDIR/$SAMPLE

    echo "[$(date)] === Starting: $SAMPLE ===" | tee -a $LOGDIR/align_all.log

    $STAR \
        --runThreadN 20 \
        --genomeDir $INDEX \
        --readFilesIn $R1 $R2 \
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
        --outFileNamePrefix $OUTDIR/$SAMPLE/ \
        >> $LOGDIR/align_all.log 2>&1

    if [ $? -eq 0 ]; then
        echo "[$(date)] === Finished successfully: $SAMPLE ===" | tee -a $LOGDIR/align_all.log
    else
        echo "[$(date)] === FAILED: $SAMPLE ===" | tee -a $LOGDIR/align_all.log
    fi

done

echo "=== All alignments complete $(date) ===" | tee -a $LOGDIR/align_all.log
