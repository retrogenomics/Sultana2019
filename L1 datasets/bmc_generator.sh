####### load preferences and define global settings
script_name="bmc_generator"
script_version='1.0'

INPUT_FILE=""
OUTPUT_FILE=""
GENOME=""

# store usage explanations
USAGE="\
$script_name v$script_version:\tStarting from a .bed file with intervals of fixed length (short),\n\
generates a random .bed file with intervals of identical size and matching base composition.\n\
usage:\t$( basename $0 ) [options] -g <genome.fa> -i <input.bed> -o <output.bed>\n\
options:\n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
"

# parse script arguments
while getopts 'hvg:i:o:' opt ; do
	case $opt in
		g) GENOME=$OPTARG ;;
		i) INPUT_FILE=$OPTARG ;;
		o) OUTPUT_FILE=$OPTARG ;;
		h) echo -e "\n$USAGE"; exit 1 ;;
		v) echo -e "${script_name} v${script_version}" ; exit 1 ;;
		\?) echo -e "\nInvalid option: -$OPTARG\n" >&2; echo -e $USAGE; exit 1 ;;
	esac
done

# check for mandatory positional parameters

if [[ -z "${INPUT_FILE}" || ! -f "${INPUT_FILE}" ]];
then
	echo -e "\nInput file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

if [[ -z "${GENOME}" || ! -f "${GENOME}" ]];
then
	echo -e "\nReference GENOME sequence (.fa) not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

# define global variables
WORKING_DIR=$( pwd )
INPUT_NAME=$( basename ${INPUT_FILE} .bed )
GENOME_NAME=$( basename ${GENOME} .fa )
GENOME_DIR=$( dirname ${GENOME} )

if [[ -z "${GENOME_DIR}/${GENOME_NAME}.genome" || ! -f "${GENOME_DIR}/${GENOME_NAME}.genome" ]];
then
	echo -e "\nReference genome size file (.genome) not specified or not existing.\n";
	echo -e "${GENOME_DIR}/${GENOME_NAME}.genome"
	echo -e $USAGE ; exit 1
fi

# calculate number of lines in input file
NB_INSERTIONS=$( awk '$1!~/^#/' "${INPUT_FILE}" | wc -l )

# calculate input file interval size and verify that it is fixed
INTERVAL_COUNT=$( awk '$1!~/^#/ {print $3-$2}' "${INPUT_FILE}" | sort -k1,1n | uniq | wc -l)
if [[ "${INTERVAL_COUNT}" -gt 1 ]];
	then
		echo -e "\nInterval length in ${INPUT_FILE} is not fixed.\n";
		echo -e $USAGE ; exit 1
	else
		l=$( awk '$1!~/^#/ {print $3-$2}' "${INPUT_FILE}" | sort -k1,1n | uniq )
		GC_WINDOW=$l
fi

# calculate %GC for each interval of the input file
echo -ne "Calculate %GC for each interval of the input file..."
if [[ -z "${INPUT_NAME}.withGCcontent.bed" || ! -f "${INPUT_NAME}.withGCcontent.bed" ]];
then
	bedtools nuc -fi ${GENOME} -bed "${INPUT_FILE}" \
	| awk '$1!~/^#/{for (i=1;i<=6;i++) {printf $i "\t"}; printf "%.02f\n",$8}' \
	| sort -k1,1 -k2,2n \
	> "${INPUT_NAME}.withGCcontent.bed"
	echo -e "Done"
else
	echo -e "Already done"
fi

# Store all the possible %GC found in input file and their number of occurence
echo -ne "Store all the possible %GC found in input file and their number of occurence..."
list_percent=$( awk '$1!~/^#/ {print $7}' "${INPUT_NAME}.withGCcontent.bed" | sort -k1,1n | uniq | awk -F"\n" '{printf $1" "}' | sed 's/.$//' )
declare -A nb
for k in ${list_percent};
do
	nb[${k}]=$( awk -v k=$k '$NF==k' "${INPUT_NAME}.withGCcontent.bed" | wc -l )
done
echo -e "Done"

# Generate random insertions matching GC content of input file
echo -ne "Generate random insertions matching GC content of input file..."
r=${NB_INSERTIONS}
tmp=""
while [ "$r" -gt 0 ] ;
do
	newline=$( bedtools random -l ${GC_WINDOW} -n 1 -g "${GENOME_DIR}/${GENOME_NAME}.genome" \
	| bedtools nuc -fi ${GENOME} -bed - \
	| tail -1 \
	| awk '($1!~/^#/ && $13==0) {for (i=1;i<=6;i++) {printf $i "\t"}; printf "%.02f\n",$8}' )
	p=$( echo "$newline" | awk '{printf $NF}' )
	if [[ "${nb[${p}]}" -gt 0 ]] && [[ -n "$p" ]]
	then
		tmp="$tmp\n$newline"
		if (( $r % 100 == 0 ))
		then
			echo -ne "|"
		fi
		(( --r ))
		(( --nb[${p}] ))
	fi
done
echo -e "...Done"

echo -e "$tmp" | awk '$1!=""' | sort -k1,1 -k2,2n > "${OUTPUT_FILE}"
exit
