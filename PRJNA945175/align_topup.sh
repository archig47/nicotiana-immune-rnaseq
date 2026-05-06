#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

TRIMDIR=/users/fyp/fyp5/project/PRJNA945175/trimmed
OUTDIR=/users/fyp/fyp5/project/PRJNA945175/aligned
LOGDIR=/users/fyp/fyp5/project/PRJNA945175/logs
INDEX=/users/fyp/fyp5/project/genome/STAR_index_softmasked
STAR=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/STAR-avx2

SAMPLES=(SRR23875083 SRR23875084 SRR23875085 SRR23875086)

echo "=== Final top-up alignment started $(date) ===" | tee -a $LOGDIR/align_topup2.log

for SAMPLE in "${SAMPLES[@]}"; do
    R1=$TRIMDIR/${SAMPLE}_1_val_1.fq.gz
    R2=$TRIMDIR/${SAMPLE}_2_val_2.fq.gz

    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "[$(date)] ERROR: Missing trimmed files for $SAMPLE — skipping" | tee -a $LOGDIR/align_topup2.log
        continue
    fi

    mkdir -p $OUTDIR/$SAMPLE

    echo "[$(date)] === Starting: $SAMPLE ===" | tee -a $LOGDIR/align_topup2.log

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
        >> $LOGDIR/align_topup2.log 2>&1

    if [ $? -eq 0 ]; then
        echo "[$(date)] === Finished successfully: $SAMPLE ===" | tee -a $LOGDIR/align_topup2.log
        rm -f $TRIMDIR/${SAMPLE}_1_val_1.fq.gz
        rm -f $TRIMDIR/${SAMPLE}_2_val_2.fq.gz
        echo "[$(date)] Deleted trimmed FASTQs for $SAMPLE" | tee -a $LOGDIR/align_topup2.log
    else
        echo "[$(date)] === FAILED: $SAMPLE ===" | tee -a $LOGDIR/align_topup2.log
    fi

done

echo "=== Final top-up complete $(date) ===" | tee -a $LOGDIR/align_topup2.log
