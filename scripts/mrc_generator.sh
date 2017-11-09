#!/usr/bin/env bash


#################################################################################
# Generate l1 matched random control (mrc) dataset. Base composition in a window
# of w bp around center of input interval is conserved
#################################################################################

#################################################################################
# Load default folders for project, picard tools, reference genome, etc
#################################################################################

CURRENT_DIR=$( pwd )
SCRIPT_PATH="`dirname \"$0\"`"
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"

# test if CONFIG file exists
configuration_file="${SCRIPT_PATH}/CONFIG"
if [[ ! -e "$configuration_file" ]];
	then
		echo -e "\nMissing configuration file in ${CURRENT_DIR}.\n";
		exit 1
fi

# read CONFIG file
while read line
do
    eval $( awk '$1!~/^#/ {print $0}' )
done < "${configuration_file}"

#################################################################################
# Set global parameters, variables and folders
#################################################################################

script_name="random_generator"
script_version='1.0'

GENOME=""
GENOME_SEQ=""
ALLOWED=""
INPUT_FILE=""
OUTPUT_DIR="${CURRENT_DIR}/mrc_dataset"
OUTPUT_FILE_PREFIX="mrc."
OUTPUT_FILE_SUFFIX=""
INTERVAL_NAME="mrc"
BOOTSTRAP=10
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
while getopts 'hvg:f:a:i:o:n:p:s:N:w:' opt ; do
	case $opt in
		g) GENOME=$OPTARG ;;
		f) GENOME_SEQ=$OPTARG ;;
		a) ALLOWED=$OPTARG ;;
		i) INPUT_FILE=$OPTARG ;;
		o) OUTPUT_DIR=$OPTARG ;;
		p) OUTPUT_FILE_PREFIX=$OPTARG ;;
		s) OUTPUT_FILE_SUFFIX=$OPTARG ;;
		n) INTERVAL_NAME=$OPTARG ;;
		N) BOOTSTRAP=$OPTARG ;;
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
# Adjust the length of interval spanning each insertion to the size of GC_WINDOW
#################################################################################

# to modify to allow INTERVAL_LENGTH<>2

echo -ne "Extract central window of ${GC_WINDOW}bp from each intervals of the input file ..."

length=$( echo ${GC_WINDOW} | awk -v i=${INTERVAL_LENGTH} '{print int(($1/2)+0.5)-(i/2)}' )
bedtools slop -b $length -i "${INPUT_FILE}" -g "${GENOME}" \
> "tmp.${GC_WINDOW}bp.${INPUT_FILE}"

echo -e "Done"

#################################################################################
# Calculate %GC for each interval of the input file if not previously done
#################################################################################

echo -ne "Calculate %GC for each adjusted intervals of the input file..."

bedtools nuc -fi "${GENOME_SEQ}" -bed "tmp.${GC_WINDOW}bp.${INPUT_FILE}" \
| awk '$1!~/^#/{for (i=1;i<=6;i++) {printf $i "\t"}; printf "%.02f\n",$8}' \
| sort -k1,1 -k2,2n \
> "tmp.withGCcontent.${INPUT_FILE}"

echo -e "Done"


# run parallel instances of mrc() function
mkdir -p ${OUTPUT_DIR}
# script_start="parallel "${SCRIPTS}/mrc_generator_single.sh" -i "tmp.withGCcontent.${INPUT_FILE}" -a ${ALLOWED} -g ${GENOME} -f ${GENOME_SEQ} -o "${OUTPUT_DIR}/{}.tmp.${GC_WINDOW}.${INPUT_FILE}" -w ${GC_WINDOW} ::: $( printf "{%04d..%04d}" 1 ${BOOTSTRAP} )"
script_start=""${SCRIPTS}/mrc_generator_single.sh" -i "tmp.withGCcontent.${INPUT_FILE}" -a ${ALLOWED} -g ${GENOME} -f ${GENOME_SEQ} -o "${OUTPUT_DIR}/0000.tmp.${GC_WINDOW}.${INPUT_FILE}" -w ${GC_WINDOW} "
eval ${script_start}

# modify coordinates of mrc to span only 2nt-intervals
for file in ${OUTPUT_DIR}/*.tmp.${GC_WINDOW}.${INPUT_FILE};
do
	file_index=$(basename "${file}" .tmp.${GC_WINDOW}.${INPUT_FILE})
	awk -v OFS="\t" -v l=${length} -v name=${INTERVAL_NAME} 'BEGIN{i=1} ($1!~/^#/) {printf $1 "\t" $2+l "\t" $2+l+2 "\t" name "|%04d\t1\t" $6 "\n", i; i++}' $file \
	> "${OUTPUT_DIR}/${OUTPUT_FILE_PREFIX}${file_index}${OUTPUT_FILE_SUFFIX}.bed"
done
