#!/opt/local/bin/bash
# update the EXP_CONFIG file and the folders below
# from a serie of atlas-seq neo run, identify insertion points and compare multiple independent experiments
# run time ~2h30

# define folders
BIOINFO="/Volumes/HD2/Lab/bioinfo"
EXP_CONFIG="${BIOINFO}/references/experiments/atlas-neo-R01-to-R09.txt"
DATA="${BIOINFO}/data/atlas-neo"
VERSION="3.2"
DATE="170418"
OUTPUT_DIR="${DATE}_atlas-neo_R01-to-R09"
CURRENT_DIR=$( pwd )
REF_GENOME_DIR="$BIOINFO/references/human"
REF_GENOME="hg19"

# make result directory
mkdir -p ${OUTPUT_DIR}

# read EXP_CONFIG file
i=1
while read aLine ;
do
	if ! [[ -z "$aLine" || "$aLine" =~ ^\# || -z "${aLine// }" ]]
	then
		run[$i]=$( echo $aLine | awk '{print $1}' )
		file[$i]=$( echo $aLine | awk '{print $2}' )
		barcode_name[$i]=$( echo $aLine | awk '{print $3}' )
		barcode_seq[$i]=$( echo $aLine | awk '{print $4}' )
		sample_name[$i]="R${run[$i]}_"$( echo $aLine | awk '{print $5}' )
		i=$(($i+1))
	fi
done < "${EXP_CONFIG}"

# identify independent sequencing runs
uniq_run=($( printf "%s\n" ${run[*]} | sort -u ))

# call atlas-clustering-neo.sh script for each run
# note that it assumes that a given run is associated to a unique .fastq file

for value in ${uniq_run[*]};
do
	# generate bc file for each run
	echo -ne "" > "${OUTPUT_DIR}/BC_R${value}.txt"
	for i in $( seq 1 ${#run[@]} );
	do
		if [[ "${run[$i]}" == "${value}" ]]; then
       		echo -e "${barcode_name[$i]}\t${barcode_seq[$i]}\t${sample_name[$i]}" \
       		>> "${OUTPUT_DIR}/BC_R${value}.txt";
       		fastq="${file[$i]}";
   		fi
   	done
  	# run identify insertion in each run
   	atlas-clustering-neo.sh -o "${OUTPUT_DIR}" -b "${OUTPUT_DIR}/BC_R${value}.txt" "${DATA}/${fastq}"
   	mv "${OUTPUT_DIR}/global.neo.3atlas.v${VERSION}.log" "${OUTPUT_DIR}/R${value}_global.neo.3atlas.v${VERSION}.log"
done

# multi-atlas-compare within each run
cd "${OUTPUT_DIR}"
for value in ${uniq_run[*]};
do
	name_list="";
	file_list="";
	for file in R${value}*.neo.3atlas.v${VERSION}.insertions.true.bed ;
	do
		name=$( echo -e $file | awk -F "." '{printf $1}' ) ;
		name_list=${name_list}$(printf ${name}",") ;
		file_list=${file_list}$(printf "${file} ") ;
	done;

	# test if a run contains a single sample and if so, small trick to run atlas compare
	nb_files=$( wc -w <<< "${file_list}" ) ;
	if [[ ${nb_files} -eq 1 ]]
	then
		name_list=${name_list}$( echo ${name_list%?} ) ;
		file_list=${file_list}$( echo ${file_list%?} ) ;
		atlas-compare-samples.sh -n ${name_list} -s NB_NRR ${file_list} ;
		cut -f1-8 multi-atlas-compare.NB_NRR.${name}.${name}.tab \
		| awk '$1~/^#/ {print $0} $1!~/^#/ {OFS="\t"; print $1,$2,$3,$4,1,$6,$7,$8}' \
		> multi-atlas-compare.NB_NRR.R${value}.tab ;
		rm multi-atlas-compare.NB_NRR.${name}.${name}.tab ;
	else
		name_list=$( echo ${name_list%?} );
		file_list=$( echo ${file_list%?} );
		atlas-compare-samples.sh -n ${name_list} -s NB_NRR ${file_list} ;
		mv multi-atlas-compare.NB_NRR.R${value}_*.tab multi-atlas-compare.NB_NRR.R${value}.tab ;
	fi

	# correct ID name and coordinates for intra-run comparison
	sed "s/EXP/R${value}/g" "multi-atlas-compare.NB_NRR.R${value}.tab" \
	| awk '{OFS="\t"; if ($1~/^#/) {print $1,$2,$3,$4,"NB_SAMPLES",$6} else {printf $1 "\t%.0f\t",($3+$2)/2; printf "%.0f\t",($3+$2)/2; printf $4 "\t" $5 "\t" $6 "\n"}}' \
	> "R${value}.insertions.true.bed" ;

	echo ;
done

# multi-atlas-compare between runs
name_list="";
file_list="";
for value in ${uniq_run[*]};
do
	name="R${value}.insertions.true.bed"
	name_list="${name_list}R${value},";
	file_list="${file_list}R${value}.insertions.true.bed ";
done;
name_list=$(echo ${name_list%?});
file_list=$( echo ${file_list%?} );
atlas-compare-samples.sh -n ${name_list} -s NB_SAMPLES ${file_list} ;

# reformat coordinates to take into account bedtools merge bug on 0-length intervals
name_list2=$( echo -ne ${name_list} | awk -F"," '{OFS="."; $NF=$NF; print $0}' )

awk '$1~/^#/ {print $0} $1!~/^#/ {printf $1"\t%.0f\t",($3+$2)/2; printf "%.0f\t",($3+$2)/2; for (i=4; i<NF; i++) {printf $i "\t"}; printf $NF "\n"}' multi-atlas-compare.NB_SAMPLES.${name_list2}.tab \
| sort -k1,1 -k2,2n \
> "R${uniq_run[0]}-R${uniq_run[-1]}.insertions.tab"

cut -f1-6 "R${uniq_run[0]}-R${uniq_run[-1]}.insertions.tab" \
| sort -k1,1 -k2,2n \
> "R${uniq_run[0]}-R${uniq_run[-1]}.insertions.true.bed"

# extract sequence flanking insertion sites with various window size
for WINDOW in 10 25 50 100 250 500 1000;
do
	bedtools slop -b ${WINDOW} -i "R${uniq_run[0]}-R${uniq_run[-1]}.insertions.true.bed" -g "${REF_GENOME_DIR}/${REF_GENOME}.genome" \
	| bedtools getfasta -s -name -fi "${REF_GENOME_DIR}/${REF_GENOME}.fa" -bed - -fo "R${uniq_run[0]}-R${uniq_run[-1]}.insertions.target.site.2x${WINDOW}.fa" ;
done

# back to initial directory
cd "${CURRENT_DIR}"

exit;
