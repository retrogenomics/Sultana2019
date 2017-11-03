#!/usr/bin/env bash

#################################################################################
# Script to call de novo l1 insertions from ATLAS-seq sequencing data
# obtained upon ectopic L1 retrotransposition (from plasmid)
#################################################################################


#################################################################################
# Set global parameters, variables and folders
#################################################################################

cmd_line="$_ ${@}" # save the options used for the log file
script_name="atlas-clustering-neo"
script_version="3.3"

start_time=`date +%s`
day=$(date +"[%d-%m-%Y] [%T]")

LC_NUMERIC_OLD=$LC_NUMERIC
export LC_NUMERIC="en_US.UTF-8"

CURRENT_DIR=$( pwd )
LINKER="GTGGCGGCCAGTATTCGTAGGAGGGCGCGTAGCATAGAACGT" # ATLAS-specific linker: RBMSL2 +T from A-tailing
starline=$( printf "*%.0s" {1..105} ) # separator for output
step=1	# store progress through the pipeline

#################################################################################
# Load default folders for project, picard tools, reference genome, etc
#################################################################################

# test if CONFIG file exists
configuration_file="${CURRENT_DIR}/CONFIG"
if [ -f "$configuration_file" ];
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
# Script argument parsing and processing
#################################################################################

# set defaults argument values
atlas5p=0
distance=10 # window (bp) for merging insertion points, 10 bp avoids to call additional insertions corresponding to sequencing errors at the breakpoint
barcode_file=""
threads=4
sampling=1
results_dir=""

# store usage explanations
USAGE="\
$script_name v${script_version}:\tanalysis of ATLAS-seq runs for ectopic L1 retrotransposition. \n\n\
usage:\t$( basename $0 ) [options] -b barcode_file.txt input_file.fastq \n\
options:\n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
\t-d Maximum distance between insertions to be merged [default=$distance] \n\
\t-s Subsampling of input fastq file (no=1; or indicate fraction of reads to consider, e.g. 0.01) [default=$sampling] \n\
\t-t Number of threads used for mapping [default=$threads]\n\
\t-o output directory [default=subdir created in current directory]\n\
\t-f To indicate a 5' ATLAS-seq experiment [not yet implemented] \n\
"

# parse script arguments
while getopts 'hvb:d:t:s:o:f' opt ; do
	case $opt in
		f) atlas5p=1 ;;
		d) distance=$OPTARG ;;
		b) barcode_file=$OPTARG ;;
		t) threads=$OPTARG ;;
		s) sampling=$OPTARG ;;
		o) results_dir=$( readlink -f "$OPTARG" ) ; mkdir -p "${results_dir}" ;;
		h) echo -e "\n$USAGE"; exit 1 ;;
		v) echo -e "${script_name} v${script_version}" ; exit 1 ;;
		\?) echo -e "\nInvalid option: -$OPTARG\n" >&2; echo -e $USAGE; exit 1 ;;
	esac
done

# skip over the processed options
shift $((OPTIND-1))
input_file="$1"

# check for mandatory positional parameters
if [[ -z "${input_file}" || ! -f "${input_file}" ]];
then
	echo -e "\nInput file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

if [[ -z "${barcode_file}" || ! -f "${barcode_file}" ]];
then
	echo -e "\nBarcode file not specified or not existing.\n";
	echo -e $USAGE ; exit 1
fi

# test if 5' or 3' ATLAS-seq experiment
if [ $atlas5p -eq 1 ];
then
	exp_name="neo.5"
	echo -e "Analysis of 5' ATLAS-seq data obtained from ectopic L1 expression has not yet been implemented."
	exit;
else
	exp_name="neo.3"
	min_amplicon_size=$(( ${#LINKER} + 25 ))
fi

# find and print the directory name for data, reference and barcode files
data_dir=$( cd "$( dirname "${input_file}" )"; pwd )
data_name=$( basename "${input_file}" )
barcode_dir=$( cd "$( dirname "${barcode_file}" )"; pwd )
barcode=$( basename "${barcode_file}" )
nb_barcodes=$( grep -v -e '^#' "${barcode_dir}/$barcode" | wc -l )

# obtain the name and sequence of barcodes and store them in bash tables
i=1
grep -v -e '^#' "${barcode_dir}/$barcode" \
> barcode_cleaned.txt

while read aLine ;
do
	barcode_name[$i]=$( echo $aLine | awk '{print $1}' )
	barcode_seq[$i]=$( echo $aLine | awk '{print $2}' )
	sample_name[$i]=$( echo $aLine | awk '{print $3}' )
	i=$(($i+1))
done < barcode_cleaned.txt
rm barcode_cleaned.txt

# create result directory with date and name of samples if not defined by user
if [[ -z "${results_dir}" || ! -d "${results_dir}" ]];
then
	results_dir="${CURRENT_DIR}/$( date +"%y%m%d_%H%M%S_" )${exp_name}atlas"
	for i in `seq 1 $nb_barcodes`
	do
		results_dir="${results_dir}_${sample_name[$i]}"
	done
	mkdir -p "${results_dir}"
fi
mkdir -p "${results_dir}/tmp"

#################################################################################
# Print header
#################################################################################

echo -ne "\n\
$starline\n\
$day ATLAS-seq ANALYSIS PIPELINE (ectopic L1 retrotransposition) v$script_version \n\
$cmd_line\n\
$starline\n\
${exp_name}' ATLAS-seq experiment \n\
Sequencing data:\t${data_dir}/${data_name} \n\
Reference genome file:\t${ref_genome_dir}/${ref_genome} \n\
Barcodes:\t\t${barcode_dir}/${barcode} \n\
Output directory:\t${results_dir} \n\
Samples: \n"

for i in $( seq 1 ${nb_barcodes} )
do
	echo -e "\t- ${barcode_name[$i]}: ${sample_name[$i]}"
done
echo -e "$starline"

# start generating statistic output file
output_stat="\
$starline\n\
$day ATLAS-seq ANALYSIS PIPELINE v${script_version} \n\
$starline\n\
${exp_name}' ATLAS-seq experiment \n\
Sequencing data:\t${data_dir}/${data_name} \n\
Reference genome file:\t${REF_GENOME_DIR}/${REF_GENOME} \n\
Barcodes:\t\t${barcode_dir}/${barcode} \n\
"

#################################################################################
# Sequencing data subsampling if requested by user
#################################################################################

printf "[Step $step - Data loading]\n"
if [ $( echo "$sampling<1" | bc -l ) -eq 1 ];
then
	input="${data_dir}/${data_name}"
	data_dir="${results_dir}/tmp"
	data_name="sub.${data_name}"
	seqtk sample "$input" "$sampling" > "${data_dir}/${data_name}"
fi
(( step++ ))

#################################################################################
# Library demultiplexing
#################################################################################

cd "${results_dir}/tmp"
printf "[Step $step - Demultiplexing]\n"

# create cutadapt adapter line
adapt=""
for i in `seq 1 $nb_barcodes`
do
	adapt+=" -g ${barcode_name[$i]}=^${barcode_seq[$i]} "
done

# split fastq file according to barcode sequence (and remove barcode) with cutadapt
cutadapt -e 0.10 -q 10 $adapt \
	--untrimmed-output=discarded_missing_bc.fastq \
	-o {name}.01.trimmed_bc.fastq \
	"${data_dir}/${data_name}" \
&> /dev/null
touch discarded_missing_bc.fastq

# output statistics
total_reads=$(( $( wc -l "${data_dir}/${data_name}" | awk '{print $1}' ) / 4 ))
discarded_missing_bc=$(( $( wc -l discarded_missing_bc.fastq | awk '{print $1}' ) / 4 ))

output_stat+="\
$starline\n\
Total processed reads:\t${total_reads} \n\
  - no barcode found:\t${discarded_missing_bc}/${total_reads} \n\
"

for i in `seq 1 $nb_barcodes`
do
	if [[ -s "${barcode_name[$i]}.01.trimmed_bc.fastq" ]];
	then
		barcode_out[$i]=$(( $( wc -l "${barcode_name[$i]}.01.trimmed_bc.fastq" | awk '{print $1}') / 4 ))
	else
		barcode_out[$i]=0
	fi
	output_stat+="  - ${sample_name[$i]} (${barcode_name[$i]}):\t${barcode_out[$i]}/${total_reads}\n"
done
output_stat+="$starline\n"

(( step++ ))

#################################################################################
# Pipeline branching depending on type of ATLAS experiment (5' or 3')
#################################################################################

if [ $atlas5p -eq 1 ];
then

	#################################################################################
	# 5' ATLAS-seq (not yet implemented)
	#################################################################################

	for i in $( seq 1 $nb_barcodes );
	do
		cd "${results_dir}/tmp"
		printf "[Step $step - Processing sample ${sample_name[$i]}]\n"
		# not (yet) implemented
		printf "Done \n"
		(( step++ ))
	done
else

	#################################################################################
	# 3' ATLAS-seq
	#################################################################################

	for i in $( seq 1 $nb_barcodes );
	do
		cd "${results_dir}/tmp"
		printf "[Step $step - Processing sample ${sample_name[$i]}]\n"

		#################################################################################
		# Trim ATLAS-seq linker at 5' end of reads (sequencing is outside to inside L1)
		#################################################################################

		printf "  - Trim ATLAS linker..."

		cutadapt -e 0.12 -q 10 -m ${min_amplicon_size} -g "^$LINKER" \
			--info-file="${barcode_name[$i]}.02b.trimmed_linker.tab" \
			--too-short-output="${barcode_name[$i]}.02c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.02d.discarded_missing_adapter.fastq" \
			-o "${barcode_name[$i]}.02a.trimmed_linker.fastq" \
			"${barcode_name[$i]}.01.trimmed_bc.fastq" \
		&> /dev/null

		# store statistics
		linker_in[$i]=${barcode_out[$i]}
		linker_out[$i]=$(( $( wc -l "${barcode_name[$i]}.02a.trimmed_linker.fastq" | awk '{print $1}') / 4 ))
		linker_too_short[$i]=$(( $( wc -l "${barcode_name[$i]}.02c.discarded_tooshort.fastq" | awk '{print $1}' ) / 4 ))
		size_selected_read_count[$i]=$(( ${linker_in[$i]} - ${linker_too_short[$i]} ))

		printf "Done \n"

		#################################################################################
		# Map trimmed ATLAS-seq reads
		#################################################################################

		printf "  - Map reads on ${REF_GENOME}..."

		# map reads with bwa-mem and select only those which are softclipped at their 3' end and which contain a stretch of A/T
		# immediately at the beginning of the clipped region
		# options: 	-t = nb of core used
		#			-M = for compatibility with Picard MarkDuplicates
		#			-v = verbosity mode
		#			-C = append FASTA/FASTQ comment to output

		bwa mem -t $threads -C -M -v 1 "${REF_GENOME_DIR}/${REF_GENOME}.fa" "${barcode_name[$i]}.02a.trimmed_linker.fastq" \
		2>/dev/null \
		| awk '($1~/^@/) || ($2==16 && $6~/^[0-9]+S/) || ($2==0 && $6~/[0-9]+S$/)' \
		| awk '	$1 ~ /^@/ {
					print $0;
				}
				$2 == 16 {
					a=match($6,/S/);
					l=substr($6,1,a-1);
					clipped=toupper(substr($10,1,l));
					if (l<=5) {patt=substr("AAAAA$",6-l)} else {patt="AAAAA$"};
					if (clipped~patt) {print $0};
				}
				$2 == 0 {
					a=match($6,/[0-9]+S$/);
					l=substr($6,a,length($6)-a);
					clipped=toupper(substr($10,length($10)-l+1));
					if (l<=5) {patt=substr("^TTTTT",1,l+1)} else {patt="^TTTTT"};
					if (clipped~patt) {print $0}
			}' \
		| samtools view -q 20 -F 260 -bu - \
		| java -Xmx8g -jar "${PATH_TO_PICARD}/picard.jar" SortSam \
			INPUT=/dev/stdin \
			OUTPUT="${barcode_name[$i]}.06.triminfo.aligned.bam" \
			SORT_ORDER=coordinate \
			QUIET=true \
			VERBOSITY=ERROR \
			VALIDATION_STRINGENCY=LENIENT \
			CREATE_INDEX=true \

		# test if the alignment .sam file is empty (only header, no reads) for samtools compatibility
		if [[ $( samtools view "${barcode_name[$i]}.06.triminfo.aligned.bam" | tail -1 ) =~ ^@.* ]];
			then
				mapped_count[$i]=0
			else
				# calculate and store the count of uniquely mapped reads after trimming of 5' linker and 3' primer-L1-pA in an array
				mapped_count[$i]=$( samtools view -c "${barcode_name[$i]}.06.triminfo.aligned.bam" )
		fi

		printf "Done \n"

		#################################################################################
		# Remove PCR duplicates
		#################################################################################

		printf "  - Remove PCR duplicates..."

		# mark and remove duplicate reads (identical start position), option java -Xmx8g related to memory usage (4g= 4Gb)
		java -Xmx8g -jar "${PATH_TO_PICARD}/picard.jar" MarkDuplicates \
			INPUT="${barcode_name[$i]}.06.triminfo.aligned.bam"\
			OUTPUT="${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" \
			METRICS_FILE="/dev/null/" \
			REMOVE_DUPLICATES=true \
			DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \
			ASSUME_SORTED=true \
			VALIDATION_STRINGENCY=LENIENT \
			QUIET=true \
			VERBOSITY=ERROR \
			CREATE_INDEX=true

		# test if the alignment .sam file after duplicate removal is not empty (only header, no reads) for samtools compatibility. The "!" after the "if" is to inverse the condition
		if ! [[ $( samtools view "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" | tail -1 ) =~ ^@.* ]];
		then
			unique_read_count[$i]=$( samtools view -c "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" )
		fi

		printf "Done \n"

		#################################################################################
		# Call L1 insertions
		#################################################################################

		printf "  - Call L1 insertions (.bed)..."

		# test if the .bam file for a given barcode exists and create bed files
		if [[ -f "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" ]];
		then
			# create a bed file with an entry per non-redundant read, INSERTION site, and sort it
			samtools view -h -q 20 -bu -F 260 "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" \
			| bedtools bamtobed -i - \
			| awk '{OFS="\t"; if ($6=="+") {print $1,$3,$3,$4,$5,"-"} else {print $1,$2,$2,$4,$5,"+"}}' \
			| sort -k1,1 -k2,2n \
			> "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed"
		else
			# create an empty .bed file if no .bam is found
			touch "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed"
		fi

		# create a header for the final bed file
		echo -e "#CHR\tINS_START\tINS_END\tINS_ID\tNB_NRR\tINS_STRAND" \
		> "${barcode_name[$i]}.11.numbered_insertions.bed"

		# merge overlapping insertion points or those distant from less than $distance and sort them
		# output the number of non-redundant read supporting each insertion and their names
		# warning: bedtools merge -s output has changed in bedtools v2.25
		# options used:
		# 	-s = force strandness (only merge reads in the same orientation)
		# 	-c = columns to operate
		# 	-o = operations to process on columns
		#	-d = max distance between 2 entries to allow merging

		if [[ -s "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed" ]];
		then
			bedtools merge -i "${barcode_name[$i]}.09a.triminfo.aligned.noduplicate.sorted.bed" -s -d $distance -c 2,6,2 -o count,distinct,collapse \
			| awk '$0!~/^#/ {\
				 printf $1 "\t" $2 "\t" $3 "\t" ;
				 printf "3ATLAS\t" ;
				 printf $5 "\t" $4 "\t" $7 "\n" ;
			}' \
			| awk '{\
				max=split($7,a,",") ;
				for (i in a) {freq[a[i]]++} ;
				printf $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" ;
				for (i=1; i<max; i++) {printf freq[a[i]] "|" a[i] ","} ;
				printf freq[a[max]] "|" a[max] "\n" ;
			}' \
			| awk -v sample=${sample_name[$i]} 'BEGIN{k=1}{\
				n=split($7,a,",") ;
				max=0 ;
				best_coord=0;
				n_best=1
				for (i=1;i<=n;i++) {
					split(a[i],b,"|") ;
					if (b[1]>max) {
						max=b[1] ;
						best_coord=b[2] ;
						n_best=1;
					}
					else if (b[1]==max && b[2]!=best_coord) {
						best_old=best_coord ;
						best_coord=(best_old+b[2]) ;
						n_best++ ;
					}
				}
				printf $1 "\t" ;
				printf "%.0f\t", best_coord/n_best ;
				printf "%.0f\t", best_coord/n_best ;
				printf sample"_3ATLAS_%.4d\t",k ;
				printf $5 "\t" $6 "\n" ;
				k++ ;
			}' \
			| sort -k1,1 -k2,2n \
			>> "${barcode_name[$i]}.11.numbered_insertions.bed"
		else
			touch "${barcode_name[$i]}.11.numbered_insertions.bed"
		fi

		printf "Done \n"

		#################################################################################
		# Calculate statistics
		#################################################################################

		# calculate and store the count of clusters
		insertion_count[$i]=$( grep -c -v -e "^#" ${barcode_name[$i]}.11.numbered_insertions.bed )

		# generate a minimal bed file for visualization purpose
		echo -e "#CHR\tINS_START\tINS_END\tINS_ID\tNB_NRR\tINS_STRAND" \
		> "${barcode_name[$i]}.12.insertions.display.bed"

		awk '$0!~/#/ {OFS="\t"; print $1,$2,$3,$4,$5,$6}' "${barcode_name[$i]}.11.numbered_insertions.bed" \
		>> "${barcode_name[$i]}.12.insertions.display.bed"

		# generates the final result files
		printf "  - Calculate statistics (.log)..."

		output_stat[$i]="\
	Sample:\t\t\t${sample_name[$i]} \n\
	Barcode:\t\t\t${barcode_name[$i]} \n\
	Reads in the run:\t\t\t${total_reads} \n\
	Barcode sorting:\t${total_reads}\t->\t${barcode_out[$i]} \n\
	Size selection:\t${barcode_out[$i]}\t->\t${size_selected_read_count[$i]} \n\
	Linker trimming:\t${size_selected_read_count[$i]}\t->\t${linker_out[$i]} \n\
	Unambigously mapped:\t${linker_out[$i]}\t->\t${mapped_count[$i]} \n\
	Non-redundant reads:\t${mapped_count[$i]}\t->\t${unique_read_count[$i]} \n\
	Insertions:\t${unique_read_count[$i]}\t->\t${insertion_count[$i]} \n\
	"
		echo -e "${output_stat[$i]}" > "${barcode_name[$i]}.log"

		printf "Done \n"

		#################################################################################
		# Generate files for visualization purposes
		#################################################################################

		# generate a bedgraph file for visualisation purpose
		printf "  - Compute coverage (.bedgraph)..."

		# plus strand bedgraph
		bedgraph_head_p="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(+) description="" visibility=full color=0,150,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_p" \
		> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph"

		bedtools genomecov -ibam "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" -bg -strand + \
		>> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph"

		# minus strand bedgraph
		bedgraph_head_m="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(-) description="" visibility=full color=0,150,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_m" \
		> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph"

		bedtools genomecov -ibam "${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" -bg -strand - \
		>> "${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph"

		printf "Done \n"

		#################################################################################
		# Obtain sequence around insertion sites with various windows
		#################################################################################

		printf "  - Extract target site flanking sequence..."

		bedtools slop -b 10 -i "${barcode_name[$i]}.12.insertions.display.bed" -g ${REF_GENOME_DIR}/${REF_GENOME}.genome \
		| bedtools getfasta -s -name -fi ${REF_GENOME_DIR}/${REF_GENOME}.fa -bed - -fo "${barcode_name[$i]}.13.target.site.2x10.fa"

		bedtools slop -b 25 -i "${barcode_name[$i]}.12.insertions.display.bed" -g ${REF_GENOME_DIR}/${REF_GENOME}.genome \
		| bedtools getfasta -s -name -fi ${REF_GENOME_DIR}/${REF_GENOME}.fa -bed - -fo "${barcode_name[$i]}.13.target.site.2x25.fa"

		bedtools slop -b 50 -i "${barcode_name[$i]}.12.insertions.display.bed" -g ${REF_GENOME_DIR}/${REF_GENOME}.genome \
		| bedtools getfasta -s -name -fi ${REF_GENOME_DIR}/${REF_GENOME}.fa -bed - -fo "${barcode_name[$i]}.13.target.site.2x50.fa"

		bedtools slop -b 100 -i "${barcode_name[$i]}.12.insertions.display.bed" -g ${REF_GENOME_DIR}/${REF_GENOME}.genome \
		| bedtools getfasta -s -name -fi ${REF_GENOME_DIR}/${REF_GENOME}.fa -bed - -fo "${barcode_name[$i]}.13.target.site.2x100.fa"

		bedtools slop -b 250 -i "${barcode_name[$i]}.12.insertions.display.bed" -g ${REF_GENOME_DIR}/${REF_GENOME}.genome \
		| bedtools getfasta -s -name -fi ${REF_GENOME_DIR}/${REF_GENOME}.fa -bed - -fo "${barcode_name[$i]}.13.target.site.2x250.fa"

		bedtools slop -b 1000 -i "${barcode_name[$i]}.12.insertions.display.bed" -g ${REF_GENOME_DIR}/${REF_GENOME}.genome \
		| bedtools getfasta -s -name -fi ${REF_GENOME_DIR}/${REF_GENOME}.fa -bed - -fo "${barcode_name[$i]}.13.target.site.2x1000.fa"

		printf "Done \n"

		#################################################################################
		# Move and rename files to keep
		#################################################################################

		printf "  - Clean up..."

		cd "${results_dir}"
		mv "tmp/${barcode_name[$i]}.log" "./${sample_name[$i]}.${exp_name}atlas.log"
		mv "tmp/${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" "./${sample_name[$i]}.${exp_name}atlas.bam"
		mv "tmp/${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bai" "./${sample_name[$i]}.${exp_name}atlas.bai"
		mv "tmp/${barcode_name[$i]}.12.insertions.display.bed" "./${sample_name[$i]}.${exp_name}atlas.insertions.true.bed"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x10.fa" "./${sample_name[$i]}.${exp_name}atlas.target.site.2x10.fa"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x25.fa" "./${sample_name[$i]}.${exp_name}atlas.target.site.2x25.fa"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x50.fa" "./${sample_name[$i]}.${exp_name}atlas.target.site.2x50.fa"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x100.fa" "./${sample_name[$i]}.${exp_name}atlas.target.site.2x100.fa"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x250.fa" "./${sample_name[$i]}.${exp_name}atlas.target.site.2x250.fa"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x1000.fa" "./${sample_name[$i]}.${exp_name}atlas.target.site.2x1000.fa"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.plus.bedgraph"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.minus.bedgraph"

		(( step++ ))

		printf "Done \n"

	done
fi

#################################################################################
# Output logs and stats
#################################################################################

# display log files
cd "${results_dir}"
for i in $( seq 1 $nb_barcodes );
do
	output_stat+="\
$( cat "${sample_name[$i]}.${exp_name}atlas.log" ) \n\
$starline \n\
"
done

# calculate runtime for the whole pipeline
end_time=`date +%s`
runtime=$( date -u -d @$(( end_time - start_time )) +"%T" )
day=$(date +"[%d-%m-%Y] [%T]")

output_stat+="\
$day \tRunning time: $runtime (hh:mm:ss) \n\
$starline \n\
"

echo -e "$output_stat" | tee "global.${exp_name}atlas.log"
cat "$0" > "$( basename "$0" ).script.log"

#################################################################################
# Restore system locale and delete temporary files
#################################################################################

export LC_NUMERIC=$LC_NUMERIC_OLD
rm -r tmp
exit;
