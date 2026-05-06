#!/bin/bash
# Compress all uncompressed fastq files in raw_fastq/
cd ~/project/GSE201377/raw_fastq
for f in *.fastq; do
    if [ -f "$f" ]; then
        echo "Compressing $f..."
        gzip "$f" &
    fi
done
wait
echo "All compression jobs done"
