#!/usr/bin/env bash

# Reformating insertions site (is) datasets and generating their respective matched random controls (mrc)
# file name convention: <reference genome>.<element>.<fragmentation>.<insertions or controls>.<cell type>.bed
# Insertions coordinates are described as 2bp-intervals containing the 2 flanking nucleotides.

root_dir="$HOME"
iss_folder="${root_dir}/projects/iss"
published="${iss_folder}/datasets/others/published"
processed="${iss_folder}/datasets/others/processed"
annotations="${iss_folder}/annotations"

mkdir -p "${processed}"

#################################################################################################
# sb (sleeping beauty)
# dataset from Gogol-Doring A et al. Mol Ther. 2016
# available here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE58744&format=file
#################################################################################################

# select specifically insertions obtained from sonicated DNA
# convert into 0-based coordinates
# and add a virtual positive strand (experiment was unstranded)

cd "${processed}"

awk -v OFS="\t" '$4~/soni/ {print $1,$2-1,$2+1,"sb|cd4|soni|ins","0","+"}' "${published}/GSM1419001_sites.SB.tab" \
| sort -k1,1 -k2,2n \
> hg18.tmp.bed

# convert into hg19 coordinates
liftOver hg18.tmp.bed "${root_dir}/references/human/hg18ToHg19.over.chain" hg19.tmp.bed unMapped.bed

# sort, filter and number insertions
sort -k1,1 -k2,2n hg19.tmp.bed \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| awk -v OFS="\t" 'BEGIN {i=1} {printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.sb.soni.ins.cd4.bed

# cleanup
rm unMapped.bed hg18.tmp.bed hg19.tmp.bed

#################################################################################################
# pb (piggybac)
# dataset from Gogol-Doring A et al. Mol Ther. 2016
# available here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE58744&format=file
#################################################################################################

# select specifically insertions obtained from sonicated DNA
# convert into 0-based coordinates
# and add a virtual strand (experiment was unstranded)

awk -v OFS="\t" '$4~/soni/ {print $1,$2-1,$2+1,"pb|cd4|soni|ins","0","+"}' "${published}/GSM1419000_sites.PB.tab" \
| sort -k1,1 -k2,2n \
> hg18.tmp.bed

# convert into hg19 coordinates
liftOver hg18.tmp.bed "${root_dir}/references/human/hg18ToHg19.over.chain" hg19.tmp.bed unMapped.bed

# sort, filter and number insertions
sort -k1,1 -k2,2n hg19.tmp.bed \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| awk -v OFS="\t" 'BEGIN {i=1} {printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.pb.soni.ins.cd4.bed

# cleanup
rm unMapped.bed hg18.tmp.bed hg19.tmp.bed

#################################################################################################
# mlv (murine leukemia virus)
# dataset from LaFave et al. NAR. 2014
# available here: https://research.nhgri.nih.gov/software/GeIST/download.shtml
#################################################################################################

# convert into 2 bp-intervals and bed6 format, filter, and number insertions
awk -v OFS="\t" '$1!~/^track/ {print $1,$2-1,$2+1,"mlv|k562|msei_mix|ins","0",$6}' "${published}/mlv_k562_LaFave_et_al.bed" \
| sort -k1,1 -k2,2n \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| sort -k1,1 -k2,2n \
| awk -v OFS="\t" 'BEGIN{i=1}{printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.mlv.msei_mix.ins.k562.bed

# same without recurrent insertions
awk -v OFS="\t" '$1!~/^track/ {print $1,$2-1,$2+1,"mlv|k562|msei_mix|loc","0",$6}' "${published}/mlv_k562_LaFave_et_al.bed" \
| sort -k1,1 -k2,2n \
| uniq \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| sort -k1,1 -k2,2n \
| awk -v OFS="\t" 'BEGIN{i=1}{printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.mlv.msei_mix.loc.k562.bed

#################################################################################################
# hiv (human immunodefiency virus)
# dataset from Wang GP et al. Genome Res. 2007
# downloaded from http://microb230.med.upenn.edu/ucsc/hiv.wig.bed
#################################################################################################

# manually separate the AvrII and MseI track from the original hiv.wig.bed
# 1bp-intervals are converted to 2bp intervals
# (the base downstream of insertion was the represented in the original wig file, irrespective of orientation, tested by BLATing several reads)

# reformat wig into 2bp-interval bed file, strand is set as "+" arbitrarily
awk -v OFS="\t" '$1!~/^track/{print $1,$2-1,$3,"hiv|jurkat|msei|ins","0","+"}' "${published}/hiv.is_jurkat.msei.hg17.wig" \
| sort -k1,1 -k2,2n \
> hg17.msei.tmp.bed

awk -v OFS="\t" '$1!~/^track/{print $1,$2,$3,"hiv|jurkat|avrii|ins","0","+"}' "${published}/hiv.is_jurkat.avrii.hg17.wig" \
| sort -k1,1 -k2,2n \
> hg17.avrii.tmp.bed

# liftover from hg17 to hg19
liftOver hg17.msei.tmp.bed "${root_dir}/references/human/hg17ToHg19.over.chain" hg19.msei.tmp.bed unMapped.msei.bed
liftOver hg17.avrii.tmp.bed "${root_dir}/references/human/hg17ToHg19.over.chain" hg19.avrii.tmp.bed unMapped.avrii.bed

# sort, filter and number insertions
sort -k1,1 -k2,2n hg19.msei.tmp.bed \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| awk -v OFS="\t" 'BEGIN {i=1} {printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.hiv.msei.ins.jurkat.bed

sort -k1,1 -k2,2n hg19.avrii.tmp.bed \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| awk -v OFS="\t" 'BEGIN {i=1} {printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.hiv.avrii.ins.jurkat.bed

rm ./*.tmp.bed unMapped.*.bed

#################################################################################################
# endogenous L1 insertions in HeLaS3 cells
# dataset from Philippe C et al. eLife. 2016
# downloaded from https://doi.org/10.7554/eLife.13926.018
#################################################################################################

# reformat 0-bp intervals into 2bp intervals and filter
awk -v OFS="\t" '$1!~/^#/ {print $1,$2-1,$2+1,"l1endo|helas3|soni|ins","0",$6}' "${published}/elife-13926-supp3-v2.helas3.bed" \
| sort -k1,1 -k2,2n \
| bedtools intersect -a - -b "${annotations}/hg19.helaAllowedGenomeSpace.bed" \
| sort -k1,1 -k2,2n \
| awk -v OFS="\t" 'BEGIN{i=1}{printf $1 "\t" $2 "\t" $3 "\t" $4 "|%04d\t" $5 "\t" $6 "\n", i; i++}' \
> hg19.l1endo.soni.ins.helas3.bed
