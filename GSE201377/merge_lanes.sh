#!/bin/bash
cd ~/project/GSE201377/raw_fastq

echo "Merging lanes..."

# Psy samples
cat SRR18901630_1.fastq.gz SRR18901631_1.fastq.gz SRR18901632_1.fastq.gz SRR18901633_1.fastq.gz > ../merged/GSM6062254_Psy_rep1_R1.fastq.gz
cat SRR18901630_2.fastq.gz SRR18901631_2.fastq.gz SRR18901632_2.fastq.gz SRR18901633_2.fastq.gz > ../merged/GSM6062254_Psy_rep1_R2.fastq.gz

cat SRR18901634_1.fastq.gz SRR18901635_1.fastq.gz SRR18901636_1.fastq.gz SRR18901637_1.fastq.gz > ../merged/GSM6062253_Psy_rep2_R1.fastq.gz
cat SRR18901634_2.fastq.gz SRR18901635_2.fastq.gz SRR18901636_2.fastq.gz SRR18901637_2.fastq.gz > ../merged/GSM6062253_Psy_rep2_R2.fastq.gz

cat SRR18901638_1.fastq.gz SRR18901639_1.fastq.gz SRR18901640_1.fastq.gz SRR18901641_1.fastq.gz > ../merged/GSM6062252_Psy_rep3_R1.fastq.gz
cat SRR18901638_2.fastq.gz SRR18901639_2.fastq.gz SRR18901640_2.fastq.gz SRR18901641_2.fastq.gz > ../merged/GSM6062252_Psy_rep3_R2.fastq.gz

cat SRR18901654_1.fastq.gz SRR18901655_1.fastq.gz SRR18901656_1.fastq.gz SRR18901657_1.fastq.gz > ../merged/GSM6062248_Psy_rep4_R1.fastq.gz
cat SRR18901654_2.fastq.gz SRR18901655_2.fastq.gz SRR18901656_2.fastq.gz SRR18901657_2.fastq.gz > ../merged/GSM6062248_Psy_rep4_R2.fastq.gz

cat SRR18901658_1.fastq.gz SRR18901659_1.fastq.gz SRR18901660_1.fastq.gz SRR18901661_1.fastq.gz > ../merged/GSM6062247_Psy_rep5_R1.fastq.gz
cat SRR18901658_2.fastq.gz SRR18901659_2.fastq.gz SRR18901660_2.fastq.gz SRR18901661_2.fastq.gz > ../merged/GSM6062247_Psy_rep5_R2.fastq.gz

cat SRR18901662_1.fastq.gz SRR18901663_1.fastq.gz SRR18901664_1.fastq.gz SRR18901665_1.fastq.gz > ../merged/GSM6062246_Psy_rep6_R1.fastq.gz
cat SRR18901662_2.fastq.gz SRR18901663_2.fastq.gz SRR18901664_2.fastq.gz SRR18901665_2.fastq.gz > ../merged/GSM6062246_Psy_rep6_R2.fastq.gz

# Pta samples
cat SRR18901642_1.fastq.gz SRR18901643_1.fastq.gz SRR18901644_1.fastq.gz SRR18901645_1.fastq.gz > ../merged/GSM6062251_Pta_rep1_R1.fastq.gz
cat SRR18901642_2.fastq.gz SRR18901643_2.fastq.gz SRR18901644_2.fastq.gz SRR18901645_2.fastq.gz > ../merged/GSM6062251_Pta_rep1_R2.fastq.gz

cat SRR18901646_1.fastq.gz SRR18901647_1.fastq.gz SRR18901648_1.fastq.gz SRR18901649_1.fastq.gz > ../merged/GSM6062250_Pta_rep2_R1.fastq.gz
cat SRR18901646_2.fastq.gz SRR18901647_2.fastq.gz SRR18901648_2.fastq.gz SRR18901649_2.fastq.gz > ../merged/GSM6062250_Pta_rep2_R2.fastq.gz

cat SRR18901650_1.fastq.gz SRR18901651_1.fastq.gz SRR18901652_1.fastq.gz SRR18901653_1.fastq.gz > ../merged/GSM6062249_Pta_rep3_R1.fastq.gz
cat SRR18901650_2.fastq.gz SRR18901651_2.fastq.gz SRR18901652_2.fastq.gz SRR18901653_2.fastq.gz > ../merged/GSM6062249_Pta_rep3_R2.fastq.gz

cat SRR18901666_1.fastq.gz SRR18901667_1.fastq.gz SRR18901668_1.fastq.gz SRR18901669_1.fastq.gz > ../merged/GSM6062245_Pta_rep4_R1.fastq.gz
cat SRR18901666_2.fastq.gz SRR18901667_2.fastq.gz SRR18901668_2.fastq.gz SRR18901669_2.fastq.gz > ../merged/GSM6062245_Pta_rep4_R2.fastq.gz

cat SRR18901670_1.fastq.gz SRR18901671_1.fastq.gz SRR18901672_1.fastq.gz SRR18901673_1.fastq.gz > ../merged/GSM6062244_Pta_rep5_R1.fastq.gz
cat SRR18901670_2.fastq.gz SRR18901671_2.fastq.gz SRR18901672_2.fastq.gz SRR18901673_2.fastq.gz > ../merged/GSM6062244_Pta_rep5_R2.fastq.gz

cat SRR18901674_1.fastq.gz SRR18901675_1.fastq.gz SRR18901676_1.fastq.gz SRR18901677_1.fastq.gz > ../merged/GSM6062243_Pta_rep6_R1.fastq.gz
cat SRR18901674_2.fastq.gz SRR18901675_2.fastq.gz SRR18901676_2.fastq.gz SRR18901677_2.fastq.gz > ../merged/GSM6062243_Pta_rep6_R2.fastq.gz

echo "Done! Merged files:"
ls -lh ../merged/
