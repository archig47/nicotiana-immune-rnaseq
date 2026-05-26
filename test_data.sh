#!/bin/bash
mkdir -p data/reads
cd data/reads
fasterq-dump --split-files SRR6676954 SRR6676955 SRR6676956 SRR6676957 SRR6676958 SRR6676959
gzip *.fastq
