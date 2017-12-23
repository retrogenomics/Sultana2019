# !/bin/bash

#############################################################
# hotspots
#############################################################

if [! -d "${OUTPUT_DIR}/tmp"]; then
mkdir "${OUTPUT_DIR}/tmp"
fi

cd "${OUTPUT_DIR}/tmp"

bin_list="0.5 1 10"

for bin in ${bin_list};
do

	# 1. Create bin bed files
# 	bin_size=$( echo "scale=0; 1000000*${bin}" | bc )
	bin_size=$( echo "${bin}*1000000/1" | bc )
	bedtools makewindows -g "${DATASETS}/hg19_filtered.genome" -w ${bin_size} \
	| sort -k1,1 -k2,2n \
	| awk -v bin=$bin '
		BEGIN {
			i=1;
		}
		{
			OFS="\t";
			printf $1 "\t" $2 "\t" $3 "\t"bin"Mb_";
			printf "%.4d\n",i;
			i++;
		}' \
	> "hg19_${bin}Mb_bins.bed"

	# 2. Compute available space in each interval by intersecting HelaS3 cover data and available space
	
	bedtools intersect -a "${DATASETS}/HeLaS3_nucleotide_coverage.bedgraph"  -b "${DATASETS}/hg19.helaAllowedGenomeSpace.bed" \
	| sort -k1,1 -k2,2n \
	| ${SCRIPTS}/prepare_for_hotspots.pl --win $bin_size \
	>"ref_cover_${bin}Mb.txt"  
	
	bedtools map -a  "hg19_${bin}Mb_bins.bed" -b "ref_cover_${bin}Mb.txt" -null 0 -c 5,6 -o sum \
	> "${bin}Mb_cover.txt"

	# 2b. Remove the huge and now useless ref_cover files.
	rm "ref_cover_${bin}Mb.txt"
	
	
	# 3. calculate total number of insertion / bin
	bedtools map -a "${bin}Mb_cover.txt" -b "${DATASETS}/hg19.l1neo.soni.ins.helas3.bed" -c 5 -o count \
	| awk '
		BEGIN {
			OFS="\t";
			print "#chr","start","end","name","gapless_size","gapless_size_cover_corrected","nb_ins";
		}
		{
			print $0;
		}' \
	> "nb_ins_per_${bin}Mb_bin.txt"

	# 4. calculate average per-base insertion rate (across gap-less genome)
	nb_insertions=$(( $( wc -l "${DATASETS}/hg19.l1neo.soni.ins.helas3.bed" | awk '{print $1}' )  ))

	awk -v n=$nb_insertions '
		BEGIN {l=0}
		($1!~/^#/) {l=l+$5}
		END {printf "%.10e\n", n/l}
	' "${bin}Mb_cover.txt" \
	> "perbase_ins_rate_${bin}Mb_bins.txt"

	# run Rscript to calculate Poisson probabilities with FDR
	Rscript --vanilla "${SCRIPTS}/hotspots_poisson.r" "${OUTPUT_DIR}/tmp" "nb_ins_per_${bin}Mb_bin.txt" "perbase_ins_rate_${bin}Mb_bins.txt" "nb_ins_per_${bin}Mb_bin_poisson.txt" "${bin}Mb_cover.txt"
done

#run Rscript to make the final circos plot
Rscript --vanilla "${SCRIPTS}/circos_plot_hotspots.r"


cd "${ISS_FOLDER}"
