#!/usr/bin/env bash

################################################################
# prepare table for plotting distribution of insertion spacing #
################################################################

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

#############################################################
# L1-L1 distance distribution
#############################################################

# number of randomisation used
n=100

# spacing of de novo insertions
sort -k1,1 -k2,2n "${iss_folder}/analysis/hg19.l1neo.soni.loc.helas3.bed" \
| bedtools spacing -i - \
| awk '$7!="." {print $7}' \
| sort -k1,1n \
> denovo_dist.txt

# spacing of random datasets
echo -ne ""> random_dist.txt;
for i in `seq 1 $n`; do \
j=$( printf "%04d" $i)
sort -k1,1 -k2,2n "${iss_folder}/datasets/l1neo/random_loc/hg19.l1neo.soni.random_loc.${j}.bed" \
| bedtools spacing -i - \
| awk '$7!="." {print $7}' \
>> random_dist.txt;
done
sort -k1,1n random_dist.txt > tmp
mv tmp random_dist.txt

# spacing of matched random control (mrc) datasets
echo -ne ""> mrc_dist.txt;
for i in `seq 1 $n`; do \
j=$( printf "%04d" $i)
sort -k1,1 -k2,2n "${iss_folder}/datasets/l1neo/mrc_loc_nonredundant/hg19.l1neo.soni.mrc_loc.${j}.bed" \
| bedtools spacing -i - \
| awk '$7!="." {print $7}' \
>> mrc_dist.txt;
done
sort -k1,1n mrc_dist.txt > tmp
mv tmp mrc_dist.txt

# generate a single table for importing into R
paste denovo_dist.txt random_dist.txt mrc_dist.txt \
> InsSpacing.tab
sed -i '1il1neo\trandom\tmrc' InsSpacing.tab

rm denovo_dist.txt random_dist.txt mrc_dist.txt

# run R script for graphing and performing statistical tests
Rscript --vanilla "${iss_folder}/scripts/InsSpacing.r"
