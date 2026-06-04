#!/bin/bash
# =============================================================================
# BLASTp: PRR proteins vs Arabidopsis thaliana UniProt Swiss-Prot database
# =============================================================================
# Purpose: assign functional labels to N. benthamiana PRR genes by finding
#          their closest characterised Arabidopsis homologue. Results used
#          to annotate gene labels in the PRR dotplot figure (Figure 8).
#
# Database: Arabidopsis thaliana UniProt Swiss-Prot (16,418 proteins)
#   - Downloaded from UniProt: https://www.uniprot.org (reviewed, A. thaliana)
#   - Built with: makeblastdb -in ath_swissprot.fasta -dbtype prot \
#                             -out ath_swissprot_db
#
# Input:  prr_proteins.faa  (1,252 PRR protein sequences from NbT2T v12)
# Output: prr_vs_ath_swissprot.txt  (tabular BLAST results, format 6)
# =============================================================================

QUERY=/users/fyp/fyp5/project/genome/prr_proteins.faa
DB=/users/fyp/fyp5/project/genome/ath_blast/ath_swissprot_db
OUT=/users/fyp/fyp5/project/genome/ath_blast/prr_vs_ath_swissprot.txt

blastp \
    -query   "$QUERY" \
    -db      "$DB" \
    -out     "$OUT" \
    -outfmt  6 \
    -evalue  1e-5 \
    -max_target_seqs 1 \
    -num_threads 20

echo "Done: $OUT"
