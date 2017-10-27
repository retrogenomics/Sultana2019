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
DATASET_NAME="R01-R09"
NB_INSERTIONS=$( awk '$1!~/^#/' "${DATASET_DIR}/${DATASET_FILE}" | wc -l )
BOOTSTRAP=1000 # number of random or control dataset to generate
GC_WINDOW=10 # window of matching base composition around insertion site

# move to working directory
cd ${WORKING_DIR}

# L1 de novo insertions
# adjust the length of interval spanning each insertion to exactly 2bp.

awk -v OFS="\t" '$1~/^#/ {print $0} $1!~/^#/ {print $1,$2-1,$3+1,$4,$5,$6}' "${DATASET_DIR}/${DATASET_FILE}" \
| sort -k1,1 -k2,2n \
> "${DATASET_NAME}.2bp.bed"

# adjust the length of interval spanning each insertion to the size of GC_WINDOW
length=$( echo ${GC_WINDOW} | awk '{print int(($1/2)+0.5)-1}' )
bedtools slop -b ${length} -i "${DATASET_NAME}.2bp.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.genome" \
> "${DATASET_NAME}.10bp.bed"

# generate random dataset
mkdir -p random
for i in $( seq 1 ${BOOTSTRAP});
do
	TAG=$( printf "${WORKING_DIR}/random/${DATASET_NAME}.random_%04d.2bp.bed" $i )
	bedtools random -l 2 -n ${NB_INSERTIONS} -g "${REF_GENOME_DIR}/${REF_GENOME}.genome" \
	| sort -k1,1 -k2,2n \
	| awk -v OFS="\t" 'BEGIN{i=1} ($1!~/^#/) {printf $1 "\t" $2 "\t" $3 "\trandom_%04d\t0\t" $6 "\n", i; i++}' \
	> ${TAG}
done

# generate base composition-matched random dataset (bmc)
# base composition is matched in a window of GC_WINDOW bp

mkdir -p base_comp_matched/
script_start="parallel "${BMC_SCRIPT}" -i "${DATASET_NAME}.10bp.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.fa" -o "${WORKING_DIR}/base_comp_matched/${DATASET_NAME}.bmc_{}.bed" ::: $( printf "{%04d..%04d}" 1 ${BOOTSTRAP} )"
eval ${script_start}

# repeat bmc script for some never-ending iterations:
# script_start="parallel "${BMC_SCRIPT}" -i "${DATASET_NAME}.10bp.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.fa" -o "${WORKING_DIR}/base_comp_matched/${DATASET_NAME}.bmc_{}.bed" ::: 0001 0004 0005 0008"
# eval ${script_start}

# to get the average GC content at the insertion sites
#awk 'BEGIN{total=0; mean=0; n=0} {total+=$7; n++} END{mean=total/n; print mean}' R01-R09.insertions.10bp.withGCcontent.bed

# modify coordinates of bmc to span only 2nt-intervals
cd base_comp_matched
for file in *.bed;
do
	name=$(basename "${file}" .bed)
	mv "${name}.bed" "${name}.${GC_window}bp.bed"
	awk -v OFS="\t" -v l=${length} 'BEGIN{i=1} ($1!~/^#/) {printf $1 "\t" $2+4 "\t" $2+6 "\tbmc_%04d\t0\t" $6 "\n", i; i++}' "${name}.${GC_window}bp.bed" \
	> "${name}.2bp.bed"
done
