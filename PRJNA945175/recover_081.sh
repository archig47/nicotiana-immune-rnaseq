#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

RAWDIR=/users/fyp/fyp5/project/PRJNA945175/raw
TRIMDIR=/users/fyp/fyp5/project/PRJNA945175/trimmed
ALIGNDIR=/users/fyp/fyp5/project/PRJNA945175/aligned
LOGDIR=/users/fyp/fyp5/project/PRJNA945175/logs
INDEX=/users/fyp/fyp5/project/genome/STAR_index_softmasked
STAR=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/STAR-avx2
TRIMGALORE=/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/trim_galore
SAMPLE=SRR23875081

echo "=== Recovery pipeline for $SAMPLE started $(date) ===" | tee -a $LOGDIR/recover_081.log

echo "[$(date)] Trimming $SAMPLE..." | tee -a $LOGDIR/recover_081.log
$TRIMGALORE --paired --quality 20 --length 36 --fastqc --cores 4 --output_dir $TRIMDIR $RAWDIR/${SAMPLE}_1.fastq.gz $RAWDIR/${SAMPLE}_2.fastq.gz >> $LOGDIR/recover_081.log 2>&1 || { echo "[$(date)] FAILED: Trimming" | tee -a $LOGDIR/recover_081.log; exit 1; }
echo "[$(date)] Trimming complete" | tee -a $LOGDIR/recover_081.log

echo "[$(date)] Aligning $SAMPLE..." | tee -a $LOGDIR/recover_081.log
mkdir -p $ALIGNDIR/$SAMPLE
$STAR --runThreadN 20 --genomeDir $INDEX --readFilesIn $TRIMDIR/${SAMPLE}_1_val_1.fq.gz $TRIMDIR/${SAMPLE}_2_val_2.fq.gz --readFilesCommand zcat --sjdbOverhang 100 --outSAMtype BAM SortedByCoordinate --outSAMattributes NH HI AS NM MD --outFilterMultimapNmax 20 --alignIntronMin 21 --alignIntronMax 1000000 --quantMode GeneCounts --outSJtype Standard --outSAMunmapped None --outReadsUnmapped None --outFileNamePrefix $ALIGNDIR/$SAMPLE/ >> $LOGDIR/recover_081.log 2>&1 || { echo "[$(date)] FAILED: Alignment" | tee -a $LOGDIR/recover_081.log; exit 1; }

echo "[$(date)] Alignment complete, cleaning up..." | tee -a $LOGDIR/recover_081.log
rm -f $RAWDIR/${SAMPLE}_1.fastq.gz $RAWDIR/${SAMPLE}_2.fastq.gz
rm -f $TRIMDIR/${SAMPLE}_1_val_1.fq.gz $TRIMDIR/${SAMPLE}_2_val_2.fq.gz
echo "=== Recovery pipeline complete $(date) ===" | tee -a $LOGDIR/recover_081.log
