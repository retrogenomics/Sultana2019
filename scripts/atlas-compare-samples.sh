#!/bin/bash

####### This is a stable version. Do not modify without changing the name

####### Dependencies:
# bedtools		2.25.0			https://github.com/arq5x/bedtools2
# gnu grep/awk

####### Release history
# 1.1 release notes (07/05/2016): stable version.
#	- corrected a bug which prevented to display correct nb of samples with insertion when -b is not used
#	- modified the way insertions are sorted in the output file (decreasing number of samples with the insertion, then decreasing number of reads per sample)
# 1.0 release notes
# 	- added option to avoid using strand information (-b)
# 0.9 release notes
#	- output start, end and strand header names as in input file (cluster or insertion)
#	- EXP_ID formatted with a fixed 4 number ID (to get better sorting)
# 0.8 release notes
# 	- added name of score in the output file name
# 0.7 release notes
# 	- implemented the -s option
# 0.6 release notes
#	- files and names to compare can be processed with options
# 0.4 release notes
# 	- annotation of the script

# known issues:
#	- if -d is set to larger than 0 (default), it happens that 2 clusters from 1 sample overlap the same transcluster, which gives a list using the map option "-o distinct"

script_name="atlas-compare-samples"
script_version='1.1'

# set defaults values
unset sample_files name 
distance=0 # window (bp) for merging clusters
score_to_display="NB_NRR"
use_of_strand="-s"

# Store usage explanations
USAGE="\
$script_name v$script_version:\tcompare ATLAS-seq clusters from several samples.\n\
Clusters should be in bed format. \n\
usage:\t$( basename $0 ) [options] A.bed B.bed ... Z.bed \n\
options:\n\
\t-d Maximum distance between reads to be merged into cluster [default=$distance]\n\
\t-n Sample names used as header (comma-separated) [default=1,2,3,...]\n\
\t-s Score to display in the sample column [default=$score_to_display] \n\
\t-b Inactivate use of strand information when comparing clusters [useful for WGA samples; default=off i.e. default=stranded clusters] \n\
\t-h Print this help menu. \n\
\t-v What version of $script_name are you using? \n\
"

# parse the options
while getopts 'hvs:d:n:b' opt ; do
	case $opt in
		d) distance=$OPTARG ;;
		s) score_to_display=$OPTARG ;;
		n) name=$OPTARG ;;
		b) use_of_strand="" ;;
		h) echo -e $USAGE; exit 1 ;;
		v) echo -e "\n$script_name, version:\t$script_version" ; exit 1 ;; 
		\?) echo -e "\nInvalid option: -$OPTARG\n" >&2; echo -e $USAGE; exit 1 ;;
	esac
done

# skip over the processed options
shift $((OPTIND-1))
i=1 
while [ "$1" ]
do
	sample_files[$i]="$1"
	shift
	i=$(($i+1))
done

# check for mandatory positional parameters
if [[ ${#sample_files[@]} -lt 2 ]];
	then
		echo -e "\nNo or only one bed file was provided for cross-checking. You need at least 2 files to compare.\n";
		echo -e $USAGE ; exit 1
fi

# generate list of sample names
nb_samples=${#sample_files[@]}
if [[ -z $name ]];
	then
		for i in `seq $nb_samples`
		do
			sample_name[$i]=$i;
		done
	else
		for i in `seq $nb_samples`
		do
			sample_name[$i]=$( echo $name | awk -F "," -v i=$i '{print $i}' ) 
		done
fi

for i in $( seq $nb_samples )
do
		echo "sample $i: name=${sample_name[$i]}"
done

# find the column corresponding to the requested score type (-s option) in each file to compare
declare -a column_for_score
declare -a nb_columns
declare -a col_start_header
declare -a col_end_header
declare -a col_strand_header

for i in `seq $nb_samples`; do
	column_for_score[$i]=$( awk -v score=${score_to_display} 'BEGIN {k=0} ($1 ~ /^#/) {for (j=1;j<=NF;++j) {if ($j==score) {k=j; exit}}} END {print k}' ${sample_files[$i]} )
	col_start_header[$i]=$( awk '$1~/^#/ {print $2}' ${sample_files[$i]} )
	col_end_header[$i]=$( awk '$1~/^#/ {print $3}' ${sample_files[$i]} )
	col_strand_header[$i]=$( awk '$1~/^#/ {print $6}' ${sample_files[$i]} )
	nb_columns[$i]=$( awk 'BEGIN {n=0} {if (NF>n){n=NF}} END {print n}' ${sample_files[$i]} )
	
	if [[ $i -gt 1 ]]
		then
			if [[ ${nb_columns[$i]} -ne ${nb_columns[$(($i - 1))]} ]]
				then
					echo -e "The number of columns in ${sample_name[$i]} and ${sample_name[$(($i - 1))]} seems to be different.";
					exit 1;
			fi

			if [[ ${column_for_score[$i]} -ne ${column_for_score[$(($i - 1))]} ]]
				then
					echo -e "The name of column #${column_for_score[$i]} in ${sample_name[$i]} and ${sample_name[$(($i - 1))]} seems to be different.";
					exit 1;
			fi
			
			if [[ ${col_start_header[$i]} != ${col_start_header[$(($i - 1))]} ]]
				then
					echo -e "The headers corresponding to start position in ${sample_name[$i]} and ${sample_name[$(($i - 1))]} seem to be different. It is possible that one file corresponds to clusters and the other to insertion points.";
					exit 1;
			fi

			if [[ ${col_end_header[$i]} != ${col_end_header[$(($i - 1))]} ]]
				then
					echo -e "The headers corresponding to end position in ${sample_name[$i]} and ${sample_name[$(($i - 1))]} seem to be different. It is possible that one file corresponds to clusters and the other to insertion points.";
					exit 1;
			fi
			
			if [[ ${col_strand_header[$i]} != ${col_strand_header[$(($i - 1))]} ]]
				then
					echo -e "The headers corresponding to strand in ${sample_name[$i]} and ${sample_name[$(($i - 1))]} seem to be different. It is possible that one file corresponds to clusters and the other to insertion points.";
					exit 1;
			fi

	fi	
done

if [[ ${column_for_score[1]} -eq 0 ]]
	then
		echo -e "No column entitled \"$score_to_display\" was found in ${sample_name[$i]}. The 5th column content of all files will be used as score and displayed in the final table (BED convention)."
	else
		echo -e "Column #${column_for_score[1]} is entitled \"$score_to_display\". Its content in all files will be used as score and displayed in the final table."
fi

# generate trans-sample clusters
# Principle:
# 	- all clusters from all samples are pooled in a single bed file and we generate trans-sample clusters by fusing those overlapping
# 	- then each cluster from each sample is compared again trans-sample clusters to find which one are present in each sample


# check if tmp folder exists and create it or overwrite
current_dir=$( pwd )
tmp_dir=".compare"
if [[ -d "${current_dir}/${tmp_dir}" ]];
	then
		rm -r ${tmp_dir} ;
fi
mkdir ${tmp_dir}

# create an empty tmp file
echo "" > ${tmp_dir}/Final.0.tmp.bed

# loop to process each sample/barcode one after the other
cat ${sample_files[@]} \
| sort -k1,1 -k2,2n -k3,3n \
> "${tmp_dir}/Final.0.tmp.bed"

head -1 "${tmp_dir}/Final.0.tmp.bed" \
> "${tmp_dir}/Final.1.tmp.bed"

awk '($1 !~ /^#/){print $0}' "${tmp_dir}/Final.0.tmp.bed" \
>> "${tmp_dir}/Final.1.tmp.bed"

# test here if stranded comparison, bedtools merge output a 4th column with strand when -s option used (which needs to be removed, since it is duplicated)
if [ "${use_of_strand}" == "-s" ] ;
	then
		bedtools merge -i "${tmp_dir}/Final.1.tmp.bed" -s -c 4,5,6 -o count,sum,distinct -d $distance \
		| awk '{OFS="\t"; print $1,$2,$3,$5,$6,$7}' \
		| sort -k1,1 -k2,2n -k3,3n \
		> "${tmp_dir}/Final.2.tmp.bed" ;
	else
		bedtools merge -i "${tmp_dir}/Final.1.tmp.bed" -c 4,5,6 -o count,sum,distinct -d $distance \
		| sort -k1,1 -k2,2n -k3,3n \
		> "${tmp_dir}/Final.2.tmp.bed" ;
fi
	
echo -e "#CHR|${col_start_header[1]}|${col_end_header[1]}|EXP_ID|NB_SAMPLES|${col_strand_header[1]}|tGBrowser" > ${tmp_dir}/Final.3.tmp.bed
awk -v k=1 '{OFS="|"; printf $1"|"$2"|"$3"|"; printf("EXP_ID_%04d|",k); print $4,$6,$1":"$2"-"$3; k++}' ${tmp_dir}/Final.2.tmp.bed >> ${tmp_dir}/Final.3.tmp.bed
		# add a unique ID for each trans-sample cluster. This ID should be unique for a given pool of samples (a run for ex. or an experiment)
		
cat ${tmp_dir}/Final.3.tmp.bed | tr "|" "\t" | sort -k1,1 -k2,2n -k3,3n > ${tmp_dir}/Final.4.tmp.bed
nb_clusters=$((`wc -l ${tmp_dir}/Final.4.tmp.bed | awk '{print $1}'`-1))

echo -e "Number of samples:\t$nb_samples"
echo -e "Number of clusters:\t$nb_clusters"

# map each sample cluster to all clusters found in the experiment (trans-sample clusters)

top="EXP_ID\tGBrowse\tCHR\t${col_start_header[1]}\t${col_end_header[1]}\t${col_strand_header[1]}\tNB_SAMPLES\t"
sort_order="-k5,5nr"

for i in `seq ${nb_samples}`
do
	echo -ne "Processing ${sample_name[$i]}..."
	top="$top${sample_name[$i]}\t"
	sort -k1,1 -k2,2n -k3,3n "${sample_files[$i]}" > "${tmp_dir}/sample_tmp.bed"
	bedtools map ${use_of_strand} -c ${column_for_score[1]} -o distinct -null "0" -a "${tmp_dir}/Final.4.tmp.bed" -b "${tmp_dir}/sample_tmp.bed"  | cut -f8- > "${tmp_dir}/Final.4.$i.tmp.bed"
		# options 	-c: Specify the column from the B file to map onto intervals in A
		# 			-o: Specify the operation that should be applied to -c
		#			-null: value if no match (default is ".", here "0")
	if [ $i -eq 1 ]
	then
		cat "${tmp_dir}/Final.4.$i.tmp.bed" > "${tmp_dir}/Final.5.$i.tmp.bed"
	else	
		paste -d "\t" "${tmp_dir}/Final.5.$(($i-1)).tmp.bed" "${tmp_dir}/Final.4.$i.tmp.bed" > "${tmp_dir}/Final.5.$i.tmp.bed"
	fi
	sort_order=$sort_order" -k$((7+$i)),$((7+$i))nr "
	echo -e "Done"
done

grep -v -e ^# "${tmp_dir}/Final.4.tmp.bed" | awk '{OFS="\t"; print $4,$7,$1,$2,$3,$6,$5}' > ${tmp_dir}/tmp.bed
echo -e $top > ${tmp_dir}/Final.6.tmp.bed
paste -d "\t" ${tmp_dir}/tmp.bed "${tmp_dir}/Final.5.$nb_samples.tmp.bed" | tr "|" "\t" >> ${tmp_dir}/Final.6.tmp.bed
head -1 ${tmp_dir}/Final.6.tmp.bed > ${tmp_dir}/Final.7.tmp.bed
tail -n +2 ${tmp_dir}/Final.6.tmp.bed | sort $sort_order -t $'\t' >> ${tmp_dir}/Final.7.tmp.bed

head -1 ${tmp_dir}/Final.7.tmp.bed \
| awk '{OFS="\t"; printf "#"$3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$6"\t"$2"\t"; for(j=8;j<=NF-1;++j) printf $j"\t"; print $NF}' \
> ${tmp_dir}/transclusters.bed

awk '$1 != "EXP_ID" {OFS="\t"; printf $3"\t"$4"\t"$5"\t"$1"\t"$7"\t"$6"\t"$2"\t"; for(j=8;j<=NF-1;++j) printf $j"\t"; print $NF}' ${tmp_dir}/Final.7.tmp.bed \
| sort $sort_order -t $'\t' \
>> ${tmp_dir}/transclusters.bed

# output file
output_file_name=$( echo $name | awk -v s=${score_to_display} 'BEGIN {FS=","; OFS="."} {printf "multi-atlas-compare."s"."; for (j=1;j<=NF;++j) {printf $j"."}; printf "tab"}' )
cp ${tmp_dir}/transclusters.bed ${current_dir}/$output_file_name
head ${current_dir}/$output_file_name

# cleanup
rm -r ${tmp_dir}