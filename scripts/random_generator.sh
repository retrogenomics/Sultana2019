#!/usr/bin/env bash

#################################################################################
# Script to generate n random control dataset where intervals have
# the same size as an input dataset. The random datasets are picked
# in a specified .bed file
#################################################################################


#################################################################################
# Set global parameters, variables and folders
#################################################################################

script_name="random_generator"
script_version='1.0'

CURRENT_DIR=$( pwd )
GENOME=""
ALLOWED=""
INPUT_FILE=""
OUTPUT_DIR="${CURRENT_DIR}/random_dataset"
OUTPUT_FILE_PREFIX="random"
OUTPUT_FILE_SUFFIX=""
INTERVAL_NAME="random"
BOOTSTRAP=10

# store usage explanations
USAGE="\
$script_name v$script_version:\tStarting from an input .bed file, generates N random .bed file with intervals of identical sizes. \n\n\
usage:\t$( basename $0 ) [options] -g <genome.bed> -i <input.bed>\n\
options:\n\
\t-o Output folder [default=${OUTPUT_DIR}]. \n\
\t-p Prefix of output files (before number) [default=${OUTPUT_FILE_PREFIX}] \n\
\t-s Prefix of output files (after number) [default=none] \n\
\t-n Name of individual intervals (a unique number will be appended) [default=${INTERVAL_NAME}]. \n\
\t-N Number of a random dataset to generate [default=${BOOTSTRAP}]. \n\
\t-a A .bed file to specify the regions in the genome, which are allowed [default=none]. \n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
"

#################################################################################
# Parse script arguments and verify parameters
#################################################################################

# parse script arguments
while getopts 'hvg:a:i:o:n:p:s:N:' opt ; do
	case $opt in
		g) GENOME=$OPTARG ;;
		a) ALLOWED=$OPTARG ;;
		i) INPUT_FILE=$OPTARG ;;
		o) OUTPUT_DIR=$OPTARG ;;
		p) OUTPUT_FILE_PREFIX=$OPTARG ;;
		s) OUTPUT_FILE_SUFFIX=$OPTARG ;;
		n) INTERVAL_NAME=$OPTARG ;;
		N) BOOTSTRAP=$OPTARG ;;
		h) echo -e "\n$USAGE"; exit 1 ;;
		v) echo -e "${script_name} v${script_version}" ; exit 1 ;;
		\?) echo -e "\nInvalid option: -$OPTARG\n" >&2; echo -e $USAGE; exit 1 ;;
	esac
done

# check for mandatory positional parameters

if [[ -z "${GENOME}" || ! -f "${GENOME}" ]];
then
	echo -e "\nReference genome size file (.genome) not specified or not existing.\n";
	echo -e "${GENOME}"
	echo -e $USAGE ; exit 1
fi

if [[ -z "${INPUT_FILE}" || ! -f "${INPUT_FILE}" ]];
then
	echo -e "\nInput file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

# test if genomic space to pick up random intervals is restricted
if [[ ! -z "${ALLOWED}" ]] ;
then
	if [[ ! -f "${ALLOWED}" ]];
	then
		echo -e "\nAllowed genomic space file ${ALLOWED} not found.\n";
		echo -e $USAGE ;
		exit 1;
	else
		INCL="-incl ${ALLOWED}"
	fi
else
	INCL=""
fi

#################################################################################
# Loop to generate BOOTSTRAP random datasets
#################################################################################

echo -ne "Generate random dataset"

# randomly pick intervals
mkdir -p ${OUTPUT_DIR}
step=$(( ${BOOTSTRAP} / 10 ))

for i in $( seq 1 ${BOOTSTRAP});
do
	TAG=$( printf "${OUTPUT_DIR}/${OUTPUT_FILE_PREFIX}%04d${OUTPUT_FILE_SUFFIX}.bed" $i )
	COMMAND=$( bedtools shuffle \
		${INCL} \
		-noOverlapping \
		-i "${INPUT_FILE}" \
		-g "${GENOME}" \
	| sort -k1,1 -k2,2n \
	| awk -v OFS="\t" -v name="${INTERVAL_NAME}" 'BEGIN{j=1} ($1!~/^#/) {printf $1 "\t" $2 "\t" $3 "\t" name "|%04d\t1\t" $6 "\n", j; j++}' \
	> "${TAG}" )
	eval ${COMMAND}

	# print progression
	if (( $i % step == 0 ))
	then
		echo -ne "."
	fi
done

echo -e "Done"
