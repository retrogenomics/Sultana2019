#!/usr/bin/env bash

#################################################################
# prepare table for plotting number of insertion per chromosome #
#################################################################

# folder configuration
# root_folder="/Volumes/HD2/Lab/bioinfo"
# iss_folder="${root_folder}/results/171113_L1neo_backup"
# output_folder="$HOME/Downloads/tmp"
# bam_folder="${root_folder}/results/project_L1_denovo_insertions/HeLaS3"

root_folder="$HOME"
iss_folder="${root_folder}/projects/iss"
output_folder="${iss_folder}/analysis"
bam_folder="${root_folder}/data/HeLaS3"

mkdir -p "${output_folder}"

# actual per-chromosome bp in hg19.helaAllowedGenomeSpace.bed
bedtools coverage \
  -a "${iss_folder}/annotations/hg19.helaMainChromosomes.bed" \
  -b "${iss_folder}/annotations/hg19.helaAllowedGenomeSpace.bed" \
| cut -f1,5 \
| sort -k1,1V \
> "${output_folder}/hg19.helaAllowedGenomeSpace.ChromSizeBp.tab"

# count of HeLaS3 WGS reads per chr
samtools view -b -h -L "${iss_folder}/annotations/hg19.helaAllowedGenomeSpace.bed" "${bam_folder}/SRR627568.bam" \
| bedtools coverage \
  -counts \
  -a "${iss_folder}/annotations/hg19.helaMainChromosomes.bed" \
  -b - \
| cut -f1,4 \
| sort -k1,1V \
> "${output_folder}/hg19.helaReadCoveragePerChr.tab"

# count of target loci per chr
bedtools coverage \
  -counts \
  -a "${iss_folder}/annotations/hg19.helaMainChromosomes.bed" \
  -b "${iss_folder}/datasets/l1neo/hg19.l1neo.soni.loc.helas3.bed" \
| cut -f1,4 \
| sort -k1,1V \
> "${output_folder}/hg19.l1neo.loc.CountPerChr.tab"

# count of insertion (insertions) per chr
bedtools coverage \
  -counts \
  -a "${iss_folder}/annotations/hg19.helaMainChromosomes.bed" \
  -b "${iss_folder}/datasets/l1neo/hg19.l1neo.soni.ins.helas3.bed" \
| cut -f1,4 \
| sort -k1,1V \
> "${output_folder}/hg19.l1neo.ins.CountPerChr.tab"

# create table for R import and plotting
sort -k1,1V "${iss_folder}/annotations/hg19.helaMainChromosomes.bed" \
| paste - \
  "${output_folder}/hg19.helaAllowedGenomeSpace.ChromSizeBp.tab" \
  "${output_folder}/hg19.helaReadCoveragePerChr.tab" \
  "${output_folder}/hg19.l1neo.loc.CountPerChr.tab" \
  "${output_folder}/hg19.l1neo.ins.CountPerChr.tab" \
| cut -f1,3,5,7,9,11 \
> "${output_folder}/chromosomal_distribution.tab"
sed -i '1ichr\tref\trefGapped\twgs\tloc\tins' chromosomal_distribution.tab

# run R script to plot data and calculate statistics
Rscript --vanilla "${iss_folder}/scripts/ChrDistScatterPlot.r"
