#!/bin/bash

####### This is a stable version. Do not modify without changing the name

####### Dependencies:
# seqtk			1.0-r75-dirty	https://github.com/lh3/seqtk
# cutadapt		1.9.2.dev0		https://github.com/marcelm/cutadapt
# bwa			0.7.10-r789		https://github.com/lh3/bwa
# picard tools	1.136			http://broadinstitute.github.io/picard/
# bedtools		2.25.0			https://github.com/arq5x/bedtools2
# gnu grep/awk

####### External .atlas.conf file should be placed in home folder
# Example:

	#	# default folder and dependency files
	#	# this file should be located in the home directory folder of the user
	#
	#	# where picard tools executable files have been installed
	#	path_to_picard="/usr/local/picard"
	#
	# 	# where MIRA assembler binaries have been installed
	# 	path_to_mira="/Applications/Biotools/mira_4.0.2/bin"
	#
	#	# where useful scripts are saved
	#	path_to_useful_scripts="$HOME/Lab/bioinfo/pipelines/useful_scripts"
	#
	#	# where the reference genome files and indexes are stored
	#	ref_genome_dir="$HOME/Lab/bioinfo/references/human"
	#	ref_genome="hg19.fa"
	#
	#	# where the results are stored (a new folder is created for each analysis in the following folder)
	#	# !!!!! It should be created in the USER home folder (write rights)
	#	results_root="$HOME/Lab/bioinfo/results"

####### Release history
# 3.2 release notes (14/04/2017)
#	- added an option to manually define the output directory
#	- called command line added to log file
# 3.1 release notes (22/07/2016)
#	- fixed a bug on the calculation of the best coordinate when several insertion points are merged to a single one
# 3.0 release notes (20/07/2016)
#	- fork of atlas-clustering_v2.2.sh to analyze de novo L1 insertions from a plasmid (only 3' and only breakpoints)
#	- what it does:
# 		- demultiplex reads based on barcode
# 		- trim linker
# 		- map on hg19 with bwa mem (soft clipping allowed)
# 		- recover ONLY reads which are softclipped AND with a polyA or polyT at the junction (only based on 5 first nt at breakpoint, or less if shorter)
# 		- call the insertion point based on softclipped position for each read
# 		- merge insertion points which are distant of <10bp, keep the most frequent among them or the mean coordinate for ties

####### Known issues
#	- 5' ATLAS for plasmid-borne retrotransposition assays not yet implemented

####### Load preferences and define global settings
cmd_line="$_ ${@}" # save the options used for the log file
script_name="atlas-clustering-neo"
script_version='3.2'

start_time=`date +%s`
day=$(date +"[%d-%m-%Y] [%T]")

LC_NUMERIC_OLD=$LC_NUMERIC
export LC_NUMERIC="en_US.UTF-8"

starline=$( printf "*%.0s" {1..105} ) # separator for output
step=1	# store progress through the pipeline

# Load location of picard tools, reference genome, and result folder from external ".atlas.conf" parameter file
configuration_file="~/.atlas.conf"
if [ -f "$configuration_file" ];
	then
		echo -e "\nMissing '.atlas.conf' file in home directory.\n";
		exit 1
fi
while read line
do
    eval $( awk '$1!~/^#/ {print $0}' )
done < ~/.atlas.conf

###### script argument parsing and processing

# set defaults argument values
atlas5p=0
distance=10 # window (bp) for merging insertion points, 10 bp avoids to call additional insertions corresponding to sequencing errors at the breakpoint
barcode_file=""
threads=4
sampling=1
results_dir=""

# store usage explanations
USAGE="\
$script_name v$script_version:\tanalysis of ATLAS-seq runs for ectopic L1 retrotransposition. \n\n\
usage:\t$( basename $0 ) [options] -b barcode_file.txt input_file.fastq \n\
options:\n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
\t-f To indicate a 5' ATLAS-seq experiment  \n\
\t-d Maximum distance between insertions to be merged [default=$distance] \n\
\t-s Subsampling of input fastq file (no=1; or indicate fraction of reads to consider, e.g. 0.01) [default=$sampling] \n\
\t-t Number of threads used for mapping [default=$threads]\n\
\t-o output directory [default=subdir created in current directory]\n\
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

# define ATLAS-seq primer and expected sequences, as well as minimal amplicon size
# ATLAS specific sequences
linker='GTGGCGGCCAGTATTCGTAGGAGGGCGCGTAGCATAGAACGT' # RBMSL2 +T from A-tailing
linker_rc='ACGTTCTATGCTACGCGCCCTCCTACGAATACTGGCCGCCAC'

polyA='TTTTTTTTTTTT' # poly(A) sequence (reverse complement)

RB3PA1='CGATACCGTAAGCCGAATTG'		# between SV40 promoter (Neo cassette) and end of L1 3'UTR (LOUXXX)
RB3PA1_rc='CAATTCGGCTTACGGTATCG'
L1_3end_rc='CGAACCCTGACGTCTTTATTATACTTTAAGTTTTAGGGTACATGTGCACATTGCC' # 3' end of L1 sequence from JM101/L1.3 Delta2 (reverse complement)

RB5PA2='TGGAAATGCAGAAATCACCG' # natural L1 here, needs to be modified
L1_5end_rc='TCTTCTGCGTCGCTCACGCTGGGAGCTGTAGACCGGAGCTGTTCCTATTCGGCCATCTTGGCTCCTCCCCC' # natural L1 here, needs to be modified

if [ $atlas5p -eq 1 ];
then
	exp_name="neo.5"
	min_amplicon=$RB5PA2$L1_5end_rc$linker_rc
	min_amplicon_size=${#min_amplicon}
	echo -e "Analysis of 5' ATLAS-seq data obtained from ectopic L1 expression has not yet been implemented."
	exit;
else
	exp_name="neo.3"
	min_amplicon=$linker
	min_amplicon_size=$(( ${#min_amplicon} + 25 ))
fi

# find and print the directory name for data, reference and barcode files

data_dir=$( cd "$( dirname "${input_file}" )"; pwd )
data_name=$( basename "${input_file}" )

barcode_dir=$( cd "$( dirname "${barcode_file}" )"; pwd )
barcode=$( basename "${barcode_file}" )
nb_barcodes=$( grep -v -e '^#' "${barcode_dir}/$barcode" | wc -l )

refL1_dir=$( cd "$( dirname "${reference_L1}" )"; pwd )
refL1_name=$( basename "${reference_L1}" )

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
	results_dir="${results_root}/$( date +"%y%m%d_%H%M%S_" )${exp_name}atlas_${script_version}"
	for i in `seq 1 $nb_barcodes`
	do
		results_dir="${results_dir}_${sample_name[$i]}"
	done
	mkdir -p "${results_dir}"
fi
mkdir -p "${results_dir}/tmp"

# print header
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
Reference genome file:\t${ref_genome_dir}/${ref_genome} \n\
Barcodes:\t\t${barcode_dir}/${barcode} \n\
"

####### test if subsampling required and if yes, generate subsampled input file
printf "[Step $step - Data loading]\n"
if [ $( echo "$sampling<1" | bc -l ) -eq 1 ];
then
	input="${data_dir}/${data_name}"
	data_dir="${results_dir}/tmp"
	data_name="sub.${data_name}"
	seqtk sample "$input" "$sampling" > "${data_dir}/${data_name}"
fi
(( step++ ))

####### barcode demultiplexing using cutadapt
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

####### pipeline branching depending on type of ATLAS experiment (5' or 3')
if [ $atlas5p -eq 1 ];
then
	####### 5' ATLAS-seq
	####### process each sample one after the other
	for i in $( seq 1 $nb_barcodes );
	do
		cd "${results_dir}/tmp"
		printf "[Step $step - Processing sample ${sample_name[$i]}]\n"

		# remove L1-specific primer sequence from the reads with cutadapt and too short
		printf "  - Trim target-specific primer..."
		cutadapt -e 0.12 -q 10 -m ${min_amplicon_size} -g "^$RB5PA2" \
			--info-file="${barcode_name[$i]}.02b.primer.trimming.tab" \
			--too-short-output="${barcode_name[$i]}.02c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.02d.discarded.missing.primer.fastq" \
			-o "${barcode_name[$i]}.02a.primer.trimming.fastq" \
			"${barcode_name[$i]}.01.trimmed_bc.fastq" \
		&> /dev/null

		# store statistics
		primer_in[$i]=${barcode_out[$i]}
		primer_out[$i]=$(( $( wc -l ${barcode_name[$i]}.02a.primer.trimming.fastq | awk '{print $1}') / 4 ))
		primer_too_short[$i]=$(( $( wc -l ${barcode_name[$i]}.02c.discarded_tooshort.fastq | awk '{print $1}') / 4 ))
		size_selected_read_count[$i]=$(( ${primer_in[$i]} - ${primer_too_short[$i]} ))

		printf "Done \n"

		# remove 5' L1 sequence from the reads with cutadapt (only for reads in which an L1-specific primer was found).
		printf "  - Trim L1 5' end..."
		cutadapt -e 0.12 -q 10 -m 25 -g $L1_5end_rc \
			--info-file="${barcode_name[$i]}.03b.5prime.trimming.tab" \
			--too-short-output="${barcode_name[$i]}.03c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.03d.discarded.missing.L1.5end.fastq" \
			-o "${barcode_name[$i]}.03a.5prime.trimming.fastq" \
			"${barcode_name[$i]}.02a.primer.trimming.fastq"  \
		&> /dev/null

		# store statistics
		L1_end_in[$i]=${primer_out[$i]}
		L1_end_out[$i]=$(( $( wc -l ${barcode_name[$i]}.03a.5prime.trimming.fastq | awk '{print $1}') / 4 ))

		printf "Done \n"

		# trim ATLAS-linker at 3' end of reads with cutadapt
			# note that reads not reaching the linker sequence are discarded since they cannot be used to eliminate PCR duplicates
		printf "  - Trim ATLAS linker..."
		cutadapt -e 0.12 -q 10 -m 25 -a $linker_rc \
			--info-file="${barcode_name[$i]}.04b.3prime.trimming.tab" \
			--too-short-output="${barcode_name[$i]}.04c.discarded_tooshort.fastq" \
			--untrimmed-output="${barcode_name[$i]}.04d.discarded.missing.adapter.fastq" \
			-o "${barcode_name[$i]}.04a.3prime.trimming.fastq" \
			"${barcode_name[$i]}.03a.5prime.trimming.fastq" \
		&> /dev/null

		# store statistics
		linker_in[$i]=${L1_end_out[$i]}
		linker_out[$i]=$(( $( wc -l ${barcode_name[$i]}.04a.3prime.trimming.fastq | awk '{print $1}') / 4 ))

		printf "Done \n"

		# Reverse complement using seqtk
			# This step is necessary to remove PCR duplicates using Picard.
			# 5' ATLAS reads reaching the linker are reversed complemented to start from the linker junction.
			# Picard will next eliminate those starting from the same position (after mapping), assuming they correspond to the same DNA fragmentation/linker ligation event.
			# Note: Reads not reaching the linker are discarded, thus lowering the coverage of the 5' junctions.
		printf "  - Reverse complement reads..."
		seqtk seq -r "${barcode_name[$i]}.04a.3prime.trimming.fastq" > "${barcode_name[$i]}.04f.3prime.trimming.rc.fastq"

		printf "Done \n"

		# map reads with bwa-mem (-t:nb of core used; -M:for compatibility with Picard MarkDuplicates; -v = verbosity mode)
		printf "  - Map reads on ${ref_genome}..."
		bwa mem -t $threads -M -v 1 "${ref_genome_dir}/${ref_genome}" "${barcode_name[$i]}.04f.3prime.trimming.rc.fastq" \
		2>/dev/null \
		| samtools view -q 20 -F 260 -bu - \
		| java -Xmx8g -jar "${path_to_picard}/picard.jar" SortSam \
			INPUT=/dev/stdin \
			OUTPUT="${barcode_name[$i]}.05.aligned.bam" \
			SORT_ORDER=coordinate \
			QUIET=true \
			VERBOSITY=ERROR \
			VALIDATION_STRINGENCY=LENIENT \
			CREATE_INDEX=true \

		# test if the alignment .sam file is empty (only header, no reads) for samtools compatibility
		if [[ $( samtools view ${barcode_name[$i]}.05.aligned.bam | tail -1 ) =~ ^@.* ]];
			then
				mapped_count[$i]=0
			else
				# calculate and store the count of uniquely mapped reads after trimming of 5' linker and 3' primer-L1-pA in an array
				mapped_count[$i]=$( samtools view -c "${barcode_name[$i]}.05.aligned.bam" )
		fi
		echo

		printf "  - Map reads on ${ref_genome}...Done \n"

		# mark and remove duplicate reads (identical start position)
		printf "  - Remove PCR duplicates..."
		java -Xmx8g -jar "${path_to_picard}/picard.jar" MarkDuplicates \
			INPUT="${barcode_name[$i]}.05.aligned.bam"\
			OUTPUT="${barcode_name[$i]}.06.aligned.noduplicate.bam" \
			METRICS_FILE="/dev/null/" \
			REMOVE_DUPLICATES=true \
			DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES\
			ASSUME_SORTED=true \
			VALIDATION_STRINGENCY=LENIENT \
			QUIET=true \
			VERBOSITY=ERROR \
			CREATE_INDEX=true

		# test if the alignment .sam file after duplicate removal is not empty (only header, no reads) for samtools compatibility. The "!" after the "if" is to inverse the condition
		if ! [[ $( samtools view "${barcode_name[$i]}.06.aligned.noduplicate.bam" | tail -1 ) =~ ^@.* ]];
		then
			unique_read_count[$i]=$( samtools view -c "${barcode_name[$i]}.06.aligned.noduplicate.bam" )
		fi

		printf "Done \n"

		# generate clusters (using bamtobed then merge, with number and name of non-redundant reads)
		printf "  - Generate clusters (.bed)..."

		# test if the .bam file for a given barcode exists and create bed files
		if [[ -f "${barcode_name[$i]}.06.aligned.noduplicate.bam" ]];
			then
				# create a bed file with an entry per read and sort it
				# options:
				#	-k1,1 sort alphabetically on the first field [using -k1 only would use from field 1 to end of line] and resolve ties with second field (-k2,2n) sorted numerically.
				# This .bed file contains only uniquely mapped, non-redundant reads, for which L1-specific primer, L1 5' end and linker were trimmed, and longer than 25nt.
				bedtools bamtobed -i "${barcode_name[$i]}.06.aligned.noduplicate.bam" \
				| sort -k1,1 -k2,2n \
				> "${barcode_name[$i]}.07.aligned.noduplicate.sorted.bed"

				# create a bed file with an entry per read and sort it (including redundant reads for rpm calculation)
				samtools view -q 20 -bu -F 260 "${barcode_name[$i]}.05.aligned.bam" \
				| bedtools bamtobed -i - \
				| sort -k1,1 -k2,2n \
				> "${barcode_name[$i]}.07.aligned.sorted.bed"
			else
				# create an empty .bed file if no .bam is found
				touch "${barcode_name[$i]}.07.aligned.noduplicate.sorted.bed"
				touch "${barcode_name[$i]}.07.aligned.sorted.bed"
		fi

		# merge overlapping reads or reads distant from less than $distance into clusters and sort them
			# warning: bedtools merge -s output has changed in bedtools v2.25
			# options used:
			# 	-s = force strandness (only merge reads in the same orientation)
			# 	-c = columns to operate
			# 	-o = operations to process on columns
			#	-d = max distance between 2 entries to allow merging
		bedtools merge -i "${barcode_name[$i]}.07.aligned.noduplicate.sorted.bed" -s -d $distance -c 4,4,6 -o distinct,count,distinct \
		| awk '$1!~/^#/ {OFS="\t"; print $1,$2,$3,$5,$6,$7}' \
		| sort -k1,1 -k2,2n \
		> "${barcode_name[$i]}.08a.clusters.sorted.bed"

		# calculate nb of reads (redundant) per cluster
			# warning: bedtools coverage input file convention has changed (-a and -b are like other tools since bedtools v2.24 )
		bedtools coverage -sorted -s -counts -a "${barcode_name[$i]}.08a.clusters.sorted.bed" -b "${barcode_name[$i]}.07.aligned.sorted.bed"\
		> "${barcode_name[$i]}.08b.clusters.sorted.with_total_reads.bed"

		# calculate and store the count of clusters
		insertion_count[$i]=$( grep -v -e "^#" ${barcode_name[$i]}.08a.clusters.sorted.bed | wc -l )

		# create a header for the final bed file
			# 	NB_READS = number of redundant reads
			#	RPM = reads assigned per millions of mapped reads (redundant)
			# 	TPM = unique tags assigned per millions of mapped non-redundant reads
		echo -e "#CHR\tSTART\tEND\tCLUSTER_ID\tNB_NRR\tCLUSTER_STRAND\tNB_READS\tRPM\tTPM" \
		> "${barcode_name[$i]}.09.numbered_clusters.bed"

		# replace the read names by a unique cluster ID in the name field ($4) and calculate RPM and TPM
		awk -v k=1 -v n=${mapped_count[$i]} -v u=${unique_read_count[$i]} '{OFS="\t"; print $1,$2,$3,k,$5,$6,$7,$7*1000000/n,$5*1000000/u; k++}' "${barcode_name[$i]}.08b.clusters.sorted.with_total_reads.bed"\
		| awk -v sample=${sample_name[$i]} '{OFS="\t"; printf $1 "\t" $2 "\t" $3 "\t"; printf sample"_5ATLAS_"; printf "%.4d\t",$4; printf $5 "\t" $6 "\t" $7 "\t"; printf "%.0f\t",$8; printf "%.0f\n",$9}'\
		>> "${barcode_name[$i]}.09.numbered_clusters.bed"

		# call insertion point (for 5' ATLAS-seq, clusters are in the same orientation as L1)
		awk 'BEGIN {printf "#CHR\tINS_START\tINS_END\tCLUSTER_ID\tNB_NRR\tINS_STRAND\tNB_READS\tRPM\tTPM\n"}
			($1!~/^#/){
				OFS="\t";
				if ($6=="+") {
					printf $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6;
					for (i=7; i<=NF; i++) {printf "\t"$i}; printf "\n"}
				else {
					printf $1"\t"$2"\t"$2"\t"$4"\t"$5"\t"$6;
					for (i=7; i<=NF; i++) {printf "\t"$i}; printf "\n"}
				}' "${barcode_name[$i]}.09.numbered_clusters.bed" \
		> "${barcode_name[$i]}.11.insertions.bed"

		# generate minimal bed files for visualization purpose
		cut -f1-6 "${barcode_name[$i]}.09.numbered_clusters.bed" \
		> "${barcode_name[$i]}.10.cluster.display.bed"

		cut -f1-6 "${barcode_name[$i]}.11.insertions.bed" \
		> "${barcode_name[$i]}.12.insertions.true.bed"

		printf "Done \n"

		# generate a bedgraph file for visualisation purpose
		printf "  - Compute coverage (.bedgraph)..."

		# plus strand bedgraph
		bedgraph_head_p="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(+) description="" visibility=full color=150,0,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_p" \
		> "${barcode_name[$i]}.06.aligned.noduplicate.plus.bedgraph"
		bedtools genomecov -ibam "${barcode_name[$i]}.06.aligned.noduplicate.bam" -bg -strand + \
		>> "${barcode_name[$i]}.06.aligned.noduplicate.plus.bedgraph"

		# minus strand bedgraph
		bedgraph_head_m="track type=bedGraph name=${sample_name[$i]}.${exp_name}atlas(-) description="" visibility=full color=150,0,0 priority=20 autoScale=off alwaysZero=on maxHeightPixels=32 graphType=bar viewLimits=0:300 yLineMark=0 yLineOnOff=on windowingFunction=mean smoothingWindow=3"
		echo -e "$bedgraph_head_m" \
		> "${barcode_name[$i]}.06.aligned.noduplicate.minus.bedgraph"
		bedtools genomecov -ibam "${barcode_name[$i]}.06.aligned.noduplicate.bam" -bg -strand - \
		>> "${barcode_name[$i]}.06.aligned.noduplicate.minus.bedgraph"

		printf "Done \n"

		# generates the final result and log files and clean up the temporary files
		printf "  - Calculate statistics (.log) and cleanup..."
		output_stat[$i]="\
Sample:\t\t\t${sample_name[$i]} \n\
Barcode:\t\t\t${barcode_name[$i]} \n\
Reads in the run:\t\t\t${total_reads} \n\
Barcode sorting:\t${total_reads}\t->\t${barcode_out[$i]} \n\
Size selection:\t${barcode_out[$i]}\t->\t${size_selected_read_count[$i]} \n\
L1-specific primer trimming:\t${size_selected_read_count[$i]}\t->\t${primer_out[$i]} \n\
L1 5' end trimming:\t${L1_end_in[$i]}\t->\t${L1_end_out[$i]} \n\
Linker trimming:\t${linker_in[$i]}\t->\t${linker_out[$i]} \n\
Unambigously mapped:\t${linker_out[$i]}\t->\t${mapped_count[$i]} \n\
Non-redundant reads:\t${mapped_count[$i]}\t->\t${unique_read_count[$i]} \n\
Clusters:\t${unique_read_count[$i]}\t->\t${insertion_count[$i]} \n\
"
		echo -e "${output_stat[$i]}" > "${barcode_name[$i]}.log"

		# move and rename files to keep
		cd "${results_dir}"
		echo -e "$output_stat" > "global.${exp_name}atlas.v${script_version}.log"
		mv "tmp/${barcode_name[$i]}.log" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.log"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.bam" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bam"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.bai" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bai"
		mv "tmp/${barcode_name[$i]}.09.numbered_clusters.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.clusters.full.bed"
		mv "tmp/${barcode_name[$i]}.10.cluster.display.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.clusters.true.bed"
		mv "tmp/${barcode_name[$i]}.11.insertions.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.full.bed"
		mv "tmp/${barcode_name[$i]}.12.insertions.true.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.true.bed"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.plus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.plus.bedgraph"
		mv "tmp/${barcode_name[$i]}.06.aligned.noduplicate.minus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.minus.bedgraph"

		printf "Done \n"
		(( step++ ))
	done
else
	####### 3' ATLAS-seq
	####### process each sample one after the other
	for i in $( seq 1 $nb_barcodes );
	do
		cd "${results_dir}/tmp"
		printf "[Step $step - Processing sample ${sample_name[$i]}]\n"

		# trim ATLAS-linker at 5' end of reads with cutadapt
		printf "  - Trim ATLAS linker..."
		cutadapt -e 0.12 -q 10 -m ${min_amplicon_size} -g "^$linker" \
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

		# map reads with bwa-mem
			# options: 	-t = nb of core used
			#			-M = for compatibility with Picard MarkDuplicates
			#			-v = verbosity mode
			#			-C = append FASTA/FASTQ comment to output
		printf "  - Map reads on ${ref_genome}..."
		bwa mem -t $threads -C -M -v 1 "${ref_genome_dir}/${ref_genome}" "${barcode_name[$i]}.02a.trimmed_linker.fastq" \
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
		| java -Xmx8g -jar "${path_to_picard}/picard.jar" SortSam \
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

		# remove PCR duplicates
		printf "  - Remove PCR duplicates..."

		# mark and remove duplicate reads (identical start position), option java -Xmx8g related to memory usage (4g= 4Gb)
		java -Xmx8g -jar "${path_to_picard}/picard.jar" MarkDuplicates \
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

		# generate clusters (using bamtobed then merge, with number and name of non-redundant reads)
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

		# calculate and store the count of clusters
		insertion_count[$i]=$( grep -c -v -e "^#" ${barcode_name[$i]}.11.numbered_insertions.bed )

		# generate a minimal bed file for visualization purpose
		echo -e "#CHR\tINS_START\tINS_END\tINS_ID\tNB_NRR\tINS_STRAND" \
		> "${barcode_name[$i]}.12.insertions.display.bed"

		awk '$0!~/#/ {OFS="\t"; print $1,$2,$3,$4,$5,$6}' "${barcode_name[$i]}.11.numbered_insertions.bed" \
		>> "${barcode_name[$i]}.12.insertions.display.bed"

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

		# obtain sequence around insertion site (multiple intervals)
		printf "  - Extract target site flanking sequence..."

		bedtools slop -b 10 -i "${barcode_name[$i]}.12.insertions.display.bed" -g $ref_genome_dir/hg19.genome \
		| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo "${barcode_name[$i]}.13.target.site.2x10.fasta"

		bedtools slop -b 25 -i "${barcode_name[$i]}.12.insertions.display.bed" -g $ref_genome_dir/hg19.genome \
		| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo "${barcode_name[$i]}.13.target.site.2x25.fasta"

		bedtools slop -b 50 -i "${barcode_name[$i]}.12.insertions.display.bed" -g $ref_genome_dir/hg19.genome \
		| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo "${barcode_name[$i]}.13.target.site.2x50.fasta"

		bedtools slop -b 100 -i "${barcode_name[$i]}.12.insertions.display.bed" -g $ref_genome_dir/hg19.genome \
		| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo "${barcode_name[$i]}.13.target.site.2x100.fasta"

		bedtools slop -b 250 -i "${barcode_name[$i]}.12.insertions.display.bed" -g $ref_genome_dir/hg19.genome \
		| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo "${barcode_name[$i]}.13.target.site.2x250.fasta"

		bedtools slop -b 1000 -i "${barcode_name[$i]}.12.insertions.display.bed" -g $ref_genome_dir/hg19.genome \
		| bedtools getfasta -s -name -fi $ref_genome_dir/$ref_genome -bed - -fo "${barcode_name[$i]}.13.target.site.2x1000.fasta"

		printf "Done \n"

		# generates the final result files
		printf "  - Calculate statistics (.log) and cleanup..."

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

		# move and rename files to keep
		cd "${results_dir}"
		mv "tmp/${barcode_name[$i]}.log" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.log"
		mv "tmp/${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bam" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bam"
		mv "tmp/${barcode_name[$i]}.07.triminfo.aligned.noduplicate.bai" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.bai"
		mv "tmp/${barcode_name[$i]}.12.insertions.display.bed" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.insertions.true.bed"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x10.fasta" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.target.site.2x10.fasta"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x25.fasta" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.target.site.2x25.fasta"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x50.fasta" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.target.site.2x50.fasta"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x100.fasta" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.target.site.2x100.fasta"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x250.fasta" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.target.site.2x250.fasta"
		mv "tmp/${barcode_name[$i]}.13.target.site.2x1000.fasta" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.target.site.2x1000.fasta"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.plus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.plus.bedgraph"
		mv "tmp/${barcode_name[$i]}.07b.triminfo.aligned.noduplicate.minus.bedgraph" "./${sample_name[$i]}.${exp_name}atlas.v${script_version}.minus.bedgraph"

		printf "Done \n"
		(( step++ ))
	done
fi

####### restore system variables
export LC_NUMERIC=$LC_NUMERIC_OLD

####### display log files
cd "${results_dir}"
for i in $( seq 1 $nb_barcodes );
do
	output_stat+="\
$( cat "${sample_name[$i]}.${exp_name}atlas.v${script_version}.log" ) \n\
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

echo -e "$output_stat" | tee "global.${exp_name}atlas.v${script_version}.log"
cat "$0" > "$( basename "$0" ).script.log"

####### delete intermediate files
rm -r tmp
