#!/usr/bin/env bash

#################################################################################
# Script to generate a matched random control dataset where intervals have
# the same size and same base composition as input dataset
#################################################################################

#################################################################################
# Set global parameters, variables and folders
#################################################################################

script_name="mrc_generator_single"
script_version='1.0'

CURRENT_DIR=$( pwd )
GENOME=""
GENOME_SEQ=""
ALLOWED=""
INPUT_FILE=""
OUTPUT_FILE=""
GC_WINDOW=10

# store usage explanations
USAGE="\
$script_name v$script_version:\tStarting from a .bed file, generates a random .bed file \n\
with intervals of identical size and matching base composition in a window \n\
surrounding the center of the input file intervals (matched random control or mrc dataset).\n\n\
usage:\t$( basename $0 ) [options] -g <genome.bed> -f <genome.fa> -i <input.bed>\n\
options:\n\
\t-o Output folder [default=${OUTPUT_DIR}]. \n\
\t-p Prefix of output files (before number) [default=${OUTPUT_FILE_PREFIX}] \n\
\t-s Prefix of output files (after number) [default=none] \n\
\t-n Name of individual intervals (a unique number will be appended) [default=${INTERVAL_NAME}]. \n\
\t-N Number of a random dataset to generate [default=${BOOTSTRAP}]. \n\
\t-a A .bed file to specify the regions in the genome, which are allowed [default=none]. \n\
\t-w Window (bp) surrounding center of input interval in which %GC is calculated [default=10]. \n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
"

#################################################################################
# Parse script arguments and verify parameters
#################################################################################

# parse script arguments
while getopts 'hvg:f:a:i:o:w:' opt ; do
	case $opt in
		g) GENOME=$OPTARG ;;
		f) GENOME_SEQ=$OPTARG ;;
		a) ALLOWED=$OPTARG ;;
		i) INPUT_FILE=$OPTARG ;;
		o) OUTPUT_FILE=$OPTARG ;;
		w) GC_WINDOW=OPTARG ;;
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

if [[ -z "${GENOME_SEQ}" || ! -f "${GENOME_SEQ}" ]];
then
	echo -e "\nReference genome sequence file (.fa) not specified or not existing.\n";
	echo -e "${GENOME}"
	echo -e $USAGE ; exit 1
fi

if [[ -z "${INPUT_FILE}" || ! -f "${INPUT_FILE}" ]];
then
	echo -e "\nInput file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

if [[ -z "${OUTPUT_FILE}" ]];
then
	echo -e "\nOutput file not specified.\n";
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

# calculate input file interval size and verify that it is fixed
INTERVAL_COUNT=$( awk '$1!~/^#/ {print $3-$2}' "${INPUT_FILE}" | sort -k1,1n | uniq | wc -l)
if [[ "${INTERVAL_COUNT}" -gt 1 ]];
	then
		echo -e "\nInterval length in ${INPUT_FILE} is not fixed.\n";
		echo -e $USAGE ; exit 1
	else
		INTERVAL_LENGTH=$( awk '$1!~/^#/ {print $3-$2}' "${INPUT_FILE}" | sort -k1,1n | uniq )
fi

#################################################################################
# Store all the possible %GC found in input file and their number of occurence
#################################################################################

echo -ne "Store all the possible %GC found in input file and their number of occurence..."

list_percent=$( \
awk '$1!~/^#/ {print $7}' "${INPUT_FILE}" \
| sort -k1,1n \
| uniq \
| awk -F"\n" '{printf $1" "}' \
| sed 's/.$//' \
)

echo -e "Done"

declare -A nb
for k in ${list_percent};
do
	nb[${k}]=$( awk -v k=$k '$NF==k' "${INPUT_FILE}" | wc -l )
	echo -e "A %GC of $k was found ${nb[${k}]} times"
done

#################################################################################
# Generate random insertions matching GC content of input file
#################################################################################

echo -ne "Generate random insertions matching GC content of input file..."

# calculate number of lines in input file
r=$( awk '$1!~/^#/' "${INPUT_FILE}" | wc -l )

# calculate number of sense strand in input file
plus=$( awk '$1!~/^#/ && $6=="+"' "${INPUT_FILE}" | wc -l )

# create arbitrary intervals of GC_WINDOW size (for sense and antisense orientations)
sense="chr1\t10\t$(( 10 + ${GC_WINDOW} ))\t.\t1\t+"
antisense="chr1\t10\t$(( 10 + ${GC_WINDOW} ))\t.\t1\t-"

tmp=""
while [[ "$r" -gt 0 ]] ;
do
	if [[ "$plus" -gt 0 ]] ;
	then
		newline=$( \
		echo -e $sense \
		| bedtools shuffle -noOverlapping -incl ${ALLOWED} -i - -g ${GENOME} \
		| head -1 \
		| bedtools nuc -fi ${GENOME_SEQ} -bed - \
		| tail -1 \
		| awk '($1!~/^#/ && $13==0) {printf $1 "\t" $2 "\t" $3 "\t" "." "\t" 1 "\t" $6 "\t" "%.02f\n", $8}' \
		)
		p=$( echo "$newline" | awk '{printf $NF}' )

		# add interval if line is non empty and if the corresponding base composition is not fully represented yet in the output
		if [[ -n "$p" ]]
			then
			if [[ "${distrib[${p}]}" -gt 0 ]]
			then
				echo -e "selected"
				tmp="$tmp\n$newline"
				(( --r ))
				(( --nb[${p}] ))
				(( --plus ))
			else
				echo -e "unselected"
			fi
		fi
	else
		newline=$( \
		echo -e $antisense \
		| bedtools shuffle -noOverlapping -incl ${ALLOWED} -i - -g ${GENOME} \
		| head -1 \
		| bedtools nuc -fi ${GENOME_SEQ} -bed - \
		| tail -1 \
		| awk '($1!~/^#/ && $13==0) {printf $1 "\t" $2 "\t" $3 "\t" "." "\t" 1 "\t" $6 "\t" "%.02f\n", $8}' \
		)
		p=$( echo "$newline" | awk '{printf $NF}' )

		# add interval if line is non empty and if the corresponding base composition is not fully represented yet in the output
		if [[ -n "$p" ]]
			then
			if [[ "${nb[${p}]}" -gt 0 ]]
			then
				echo -e $newline
				tmp="$tmp\n$newline"
				(( --r ))
				(( --nb[${p}] ))
			fi
		fi
	fi
done

# reformat and sort output file before saving it
echo -e "$tmp" \
| awk '$1!=""' \
| sort -k1,1 -k2,2n \
> "${OUTPUT_FILE}"

echo -e "...Done"

exit
