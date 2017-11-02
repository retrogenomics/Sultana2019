#!/opt/local/bin/bash
# update the EXP_CONFIG file and the folders below
# from a serie of atlas-seq neo run, identify insertion points and compare multiple independent experiments
# run time ~2h30

# define folders
BIOINFO="/Volumes/HD2/Lab/bioinfo"
EXP_CONFIG="${BIOINFO}/references/experiments/atlas-neo-R01-to-R09.txt"
DATA="${BIOINFO}/data/atlas-neo"
VERSION="3.2"
DATE="171031"
OUTPUT_DIR="${DATE}_atlas-neo_R01-to-R09"
CURRENT_DIR=$( pwd )
REF_GENOME_DIR="${BIOINFO}/references/human"
REF_GENOME="hg19"
DUKE_FILTER="${BIOINFO}/annotations/hg19/other/hg19.wgEncodeDukeMapabilityRegionsExcludable.bed"
ENCODE_FILTER="${BIOINFO}/annotations/hg19/other/hg19.wgEncodeDacMapabilityConsensusExcludable.bed"
GAPLESS_GENOME="${BIOINFO}/annotations/hg19/other/hg19.helaGapLessGenome.bed"

s="************"
# make result directory
mkdir -p ${OUTPUT_DIR}

# read EXP_CONFIG file
echo -e "\n$s Read experiment metadata $s"

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

echo -e "\n$s L1 insertion calling $s"

for value in ${uniq_run[*]};
do
	echo -ne "Run ${value}..."
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
#    	atlas-clustering-neo.sh -o "${OUTPUT_DIR}" -b "${OUTPUT_DIR}/BC_R${value}.txt" "${DATA}/${fastq}"
#    	mv "${OUTPUT_DIR}/global.neo.3atlas.v${VERSION}.log" "${OUTPUT_DIR}/R${value}_global.neo.3atlas.v${VERSION}.log"
	rm "${OUTPUT_DIR}/BC_R${value}.txt"
	echo -e "Done"
done

# multi-atlas-compare within each run

echo -e "\n$s Intra-run comparisons $s"

cd "${OUTPUT_DIR}"
echo -ne "" > Blacklist_atlas.bed

for value in ${uniq_run[*]};
do
	echo -e "\n$s Run #${value}"
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
		> multi-atlas-compare.NB_NRR.R${value}.tmp ;
		rm multi-atlas-compare.NB_NRR.${name}.${name}.tab ;
	else
		name_list=$( echo ${name_list%?} );
		file_list=$( echo ${file_list%?} );
		atlas-compare-samples.sh -n ${name_list} -s NB_NRR ${file_list} ;
		mv multi-atlas-compare.NB_NRR.R${value}_*.tab multi-atlas-compare.NB_NRR.R${value}.tmp ;
	fi

	# filter for Hela main chromosomes only, exclude black listed genomic regions
	# remove insertions corresponding to mapping artefact (insertion less than 500 bp from another insertion, supported with less reads, same orientation and same run)
	n=$( head -1 multi-atlas-compare.NB_NRR.R${value}.tmp | awk '{print NF}' )

	head -1 multi-atlas-compare.NB_NRR.R${value}.tmp \
	> multi-atlas-compare.NB_NRR.R${value}.tab

	field=$( echo | awk -v n=$n '{for (i=1;i<n;i++) {printf i","}; printf n}' )

	sort -k1,1 -k2,2n multi-atlas-compare.NB_NRR.R${value}.tmp \
	| bedtools merge -s -d 500 -i - -c $field -o collapse \
	| awk -v OFS="\t" -v n=$n -v h=$h '\
		BEGIN {samples=n-7}
		($5!~/,/) {for (i=5;i<n+4;i++) {printf $i "\t"}; printf $(n+4) "\n" }
		($5~/,/) {
			best=1
			highest=0
			for (i=NF-samples+1; i<=NF; i++) {
				k=split($i,cov,",")
				for (j=1;j<=k;j++) {
					if (cov[j]>highest){highest=cov[j]; best=j}
				}
			}
			for (i=1;i<=n;i++) {
				split($(i+4),f,",")
				if (i<n) {printf f[best] "\t"} else {printf f[best] "\n"}
			}
		}
 	' \
	>> multi-atlas-compare.NB_NRR.R${value}.tab

	# create a black list of excluded insertions for inter-run comparison
	bedtools intersect -v -wa -a multi-atlas-compare.NB_NRR.R${value}.tmp -b multi-atlas-compare.NB_NRR.R${value}.tab \
	| sort -k1,1 -k2,2n \
	| cut -f1-6 \
	>> Blacklist_atlas.bed

	rm multi-atlas-compare.NB_NRR.R${value}.tmp

	# correct ID name and coordinates for intra-run comparison
	sed "s/EXP/R${value}/g" "multi-atlas-compare.NB_NRR.R${value}.tab" \
	| awk '{OFS="\t"; if ($1~/^#/) {print $1,$2,$3,$4,"NB_SAMPLES",$6} else {printf $1 "\t%.0f\t",($3+$2)/2; printf "%.0f\t",($3+$2)/2; printf $4 "\t" $5 "\t" $6 "\n"}}' \
	> "R${value}.insertions.true.bed" ;
done

# multi-atlas-compare between runs
echo -e "\n$s Inter-run comparison and insertion filtering $s"

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
| bedtools intersect -v -a - -b ${DUKE_FILTER} \
| bedtools intersect -v -a - -b ${ENCODE_FILTER} \
| bedtools intersect -v -a - -b Blacklist_atlas.bed \
| bedtools intersect -a - -b ${GAPLESS_GENOME} \
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
