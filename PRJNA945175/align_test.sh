#!/bin/bash
source /users/fyp/fyp5/miniconda3/etc/profile.d/conda.sh
conda activate nlr_pipeline

mkdir -p /users/fyp/fyp5/project/PRJNA945175/aligned/SRR23875051

echo "=== Test alignment started $(date) ===" 

/users/fyp/fyp5/miniconda3/envs/nlr_pipeline/bin/STAR-avx2 \
  --runThreadN 20 \
  --genomeDir /users/fyp/fyp5/project/genome/STAR_index_softmasked \
  --readFilesIn \
    /users/fyp/fyp5/project/PRJNA945175/trimmed/SRR23875051_1_val_1.fq.gz \
    /users/fyp/fyp5/project/PRJNA945175/trimmed/SRR23875051_2_val_2.fq.gz \
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
  --outFileNamePrefix /users/fyp/fyp5/project/PRJNA945175/aligned/SRR23875051/

echo "=== Test alignment finished $(date) ==="
