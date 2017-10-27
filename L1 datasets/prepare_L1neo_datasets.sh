#!/bin/bash

##### insertion datasets used to generate final figures of Sultana et al. article on L1 insertion site preference

# define folders
BIOINFO="$HOME"
WORKING_DIR=$( pwd )
REF_GENOME_DIR="$BIOINFO/references/human"
REF_GENOME="hg19"
BMC_SCRIPT="${WORKING_DIR}/bmc_generator.sh"

# set global variables
DATASET_FILE="R01-R09.insertions.true.bed"
DATASET_DIR="${WORKING_DIR}"
DATASET_NAME="l1neo.ins_helas3.soni.hg19"
BOOTSTRAP=1000 # number of random or control dataset to generate
GC_WINDOW=10 # window of matching base composition around insertion site

# move to working directory
cd ${WORKING_DIR}

# L1 de novo insertions
# adjust the length of interval spanning each insertion to exactly 2bp
# expand the lines representing multiple insertions at the same position

echo -ne "Reformat data files..."

sort -k1,1 -k2,2n "${DATASET_DIR}/${DATASET_FILE}" \
| awk '($1!~/^#/) {while ($5--) {print $0}}' \
| awk 'BEGIN {i=1} {printf $1 "\t" $2-1 "\t" $3+1 "\tl1neo|helas3|soni|ins|%04d\t1\t" $6 "\n", i; i++}' \
> "${DATASET_NAME}.bed"

NB_INSERTIONS=$( awk '$1!~/^#/' "${DATASET_DIR}/${DATASET_NAME}.bed" | wc -l )


# idem without expanding lines (L1 target loci)
sort -k1,1 -k2,2n "${DATASET_DIR}/${DATASET_FILE}" \
| awk 'BEGIN {i=1} ($1!~/^#/) {printf $1 "\t" $2-1 "\t" $3+1 "\tl1neo|helas3|soni|loc|%04d\t1\t" $6 "\n", i; i++}' \
> "l1neo.loc_helas3.soni.hg19.bed"

# adjust the length of interval spanning each insertion to the size of GC_WINDOW
length=$( echo ${GC_WINDOW} | awk '{print int(($1/2)+0.5)-1}' )
bedtools slop -b ${length} -i "${DATASET_NAME}.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.genome" \
> "${DATASET_NAME}.10bp.bed"

echo -e "Done"

# generate random dataset

echo -ne "Generate random dataset..."

mkdir -p random
for i in $( seq 1 ${BOOTSTRAP});
do
	TAG=$( printf "${WORKING_DIR}/random/l1neo.random_%04d.soni.hg19.bed" $i )
	bedtools random -l 2 -n ${NB_INSERTIONS} -g "${REF_GENOME_DIR}/${REF_GENOME}.genome" \
	| sort -k1,1 -k2,2n \
	| awk -v OFS="\t" 'BEGIN{i=1} ($1!~/^#/) {printf $1 "\t" $2 "\t" $3 "\tl1neo|insilico|soni|random|%04d\t1\t" $6 "\n", i; i++}' \
	> ${TAG}
done

echo -e "Done"

# generate base composition-matched random dataset (bmc)
# base composition is matched in a window of GC_WINDOW bp

mkdir -p base_comp_matched/
script_start="parallel "${BMC_SCRIPT}" -i "${DATASET_NAME}.10bp.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.fa" -o "${WORKING_DIR}/base_comp_matched/l1neo.mrc_{}.soni.hg19.${GC_WINDOW}bp.bed" ::: $( printf "{%04d..%04d}" 1 ${BOOTSTRAP} )"
eval ${script_start}

# repeat bmc script for some never-ending iterations:
# script_start="parallel "${BMC_SCRIPT}" -i "${DATASET_NAME}.10bp.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.fa" -o "${WORKING_DIR}/base_comp_matched/l1neo.bmc_{}.soni.hg19.bed" ::: 0001 0004 0005 0008"
# eval ${script_start}

# to get the average GC content at the insertion sites
#awk 'BEGIN{total=0; mean=0; n=0} {total+=$7; n++} END{mean=total/n; print mean}' l1neo.is_helas3.soni.hg19.10bp.withGCcontent.bed

# modify coordinates of bmc to span only 2nt-intervals
cd base_comp_matched
for file in *.bed;
do
	name=$(basename "${file}" .${GC_WINDOW}bp.bed)
	awk -v OFS="\t" -v l=${length} 'BEGIN{i=1} ($1!~/^#/) {printf $1 "\t" $2+l-1 "\t" $2+l+1 "\tl1neo|insilico|soni|mrc|%04d\t1\t" $6 "\n", i; i++}' $file \
	> "${name}.bed"
done
