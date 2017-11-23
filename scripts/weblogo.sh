#!/usr/bin/env bash

#################################################################################
# procedure to generate weblogo of insertion site motifs
#################################################################################

bedtools slop -s -l 2 -r 7 -i hg19.l1neo.soni.loc.helas3.bed -g ../../annotations/hg19.genome \
| bedtools getfasta -s -name -fi ~/references/human/hg19.fa -bed - -fo hg19.l1neo.soni.loc.helas3.-3bp+8bp.fa

weblogo -f hg19.l1neo.soni.loc.helas3.-3bp+8bp.fa \
-A dna \
-F pdf \
-P "" \
-c classic \
--xlabel "Position relative to insertion point" \
--annotate -3,-2,-1,1,2,3,4,5,6,7,8 \
--number-fontsize 6 \
--fontsize 8 \
--logo-font Helvetica-Bold \
--text-font Helvetica \
--aspect-ratio 5 \
-o hg19.l1neo.soni.loc.helas3.-3bp+8bp.weblogo.pdf
