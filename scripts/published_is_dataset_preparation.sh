# Reformating insertions site (is) datasets and generating their respective matched random controls (mrc)
# file name convention: <element>.<sample>.<fragmentation>.<reference genome>.bed
# where <sample> can be either insertion sites in a given cell type (e.g. is_k562) or matched random control (e.g. mrc)
# Insertions coordinates are described as 2bp-intervals containing the 2 flanking nucleotides.

#################################################################################################
# sb (sleeping beauty)
# dataset from Gogol-Döring A et al. Mol Ther. 2016
# available here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE58744&format=file
#################################################################################################

# select specifically insertions obtained from sonicated DNA
# convert into 0-based coordinates
# and add a virtual strand (experiment was unstranded)

awk -v OFS="\t" '$4~/soni/ {print $1,$2-1,$2+1,"sb|is|cd4|soni","0","+"}' GSM1419001_sites.SB.tab.txt \
| sort -k1,1 -k2,2n \
> sb.is_cd4.soni.hg18.bed

# convert into hg19 coordinates

liftOver sb.is_cd4.soni.hg18.bed ~/references/human/hg18ToHg19.over.chain sb.is_cd4.soni.hg19.bed unMapped.bed
rm unMapped.bed

# generate matched random control (any position with TA), 10 times more mrc than is

current_n=0
target_n=$(( $( wc -l sb.is_cd4.soni.hg19.bed | awk '{print $1}' ) * 10 ))
echo -ne "" > echo -ne "" >
while [ ${current_n} -ne ${target_n} ];
do
	bedtools random -l 2 -n $(( ${target_n} * 20 )) -g ~/references/human/hg19.genome \
	| bedtools getfasta -fi ~/references/human/hg19.fa -bed - -tab -s \
	| awk -F "[:()\t]" '{OFS="\t"; split($2,a,"-") ;print $1,a[1],a[2],"sb|mrc|insilico|soni","0",$3,toupper($NF)}' \
	| awk '$NF=="TA"' \
	>> tmp.bed
	current_n=$( wc -l tmp.bed | awk '{print $1}' )
	if [[ ${current_n} -ge ${target_n} ]];
	then
		head -n ${target_n} tmp.bed \
		| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6}' \
		| sort -k1,1 -k2,2n \
		> sb.mrc.soni.hg19.bed
		current_n=${target_n}
	fi
done
rm tmp.bed

#################################################################################################
# pb (piggybac)
# dataset from Gogol-Döring A et al. Mol Ther. 2016
# available here: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE58744&format=file
#################################################################################################

# select specifically insertions obtained from sonicated DNA
# convert into 0-based coordinates
# and add a virtual strand (experiment was unstranded)

awk -v OFS="\t" '$4~/soni/ {print $1,$2-1,$2+1,"pb|is|cd4|soni","0","+"}' GSM1419000_sites.PB.tab.txt \
| sort -k1,1 -k2,2n \
> pb.is_cd4.soni.hg18.bed

# convert into hg19 coordinates

liftOver pb.is_cd4.soni.hg18.bed ~/references/human/hg18ToHg19.over.chain pb.is_cd4.soni.hg19.bed unMapped.bed
rm unMapped.bed

# generate matched random control (any position with TTAA), 10 times more mrc than is

current_n=0
target_n=$(( $( wc -l pb.is_cd4.soni.hg19.bed | awk '{print $1}' ) * 10 ))
echo -ne "" > echo -ne "" >
while [ ${current_n} -ne ${target_n} ];
do
	bedtools random -l 4 -n $(( ${target_n} * 400 )) -g ~/references/human/hg19.genome \
	| bedtools getfasta -fi ~/references/human/hg19.fa -bed - -tab -s \
	| awk -F "[:()\t]" '{OFS="\t"; split($2,a,"-") ;print $1,a[1],a[2],"pb|mrc|insilico|soni","0",$3,toupper($NF)}' \
	| awk '$NF=="TTAA"' \
	>> tmp.bed
	current_n=$( wc -l tmp.bed | awk '{print $1}' )
	if [[ ${current_n} -ge ${target_n} ]];
	then
		head -n ${target_n} tmp.bed \
		| awk -v OFS="\t" '{print $1,$2+1,$3-1,$4,$5,$6}' \
		| sort -k1,1 -k2,2n \
		> pb.mrc.soni.hg19.bed
		current_n=${target_n}
	fi
done
rm tmp.bed

#################################################################################################
# mlv (murine leukemia virus)
# dataset from LaFave et al. NAR. 2014
# available here: https://research.nhgri.nih.gov/software/GeIST/download.shtml
#################################################################################################

# convert into 2 bp-intervals and bed6 format
awk -v OFS="\t" '$1!~/^track/ {print $1,$2-1,$2+1,"mlv|is|k562|msei_mix","0",$6}' mlv_k562_LaFave_et_al.bed \
| sort -k1,1 -k2,2n \
> mlv.is_k562.msei_mix.hg19.bed

#	- generate matched random control (position relative to restriction site)
random_site_driver.sh experimental_distance_from_restriction restriction_sites total_events_needed hg19.fa chr_size_UCSC_hg19 hg19_noY_noUn mlv_k562_mrc 1000

#################################################################################################
# hiv (human immunodefiency virus)
# dataset from Wang GP et al. Genome Res. 2007
# downloaded from http://microb230.med.upenn.edu/ucsc/hiv.wig.bed
#################################################################################################

# manually separate the AvrII and MseI track from the original hiv.wig.bed

# reformat wig into bed file
awk -v OFS="\t" '$1!~/^track/{print $1,$2,$3,"hiv|is|jurkat|msei","0","+"}' hiv.is_jurkat.msei.hg17.wig.bed \
| sort -k1,1 -k2,2n \
> hiv.is_jurkat.msei.hg17.bed

awk -v OFS="\t" '$1!~/^track/{print $1,$2,$3,"hiv|is|jurkat|avrii","0","+"}' hiv.is_jurkat.avrii.hg17.wig.bed \
| sort -k1,1 -k2,2n \
> hiv.is_jurkat.avrii.hg17.bed

# liftover from hg17 to hg19
liftOver hiv.is_jurkat.msei.hg17.bed ~/references/human/hg17ToHg19.over.chain tmp.bed unMapped_msei.bed
sort -k1,1 -k2,2n tmp.bed \
> hiv.is_jurkat.msei.hg19.bed

liftOver hiv.is_jurkat.avrii.hg17.bed ~/references/human/hg17ToHg19.over.chain tmp.bed unMapped_avrii.bed
sort -k1,1 -k2,2n tmp.bed \
> hiv.is_jurkat.avrii.hg19.bed

rm tmp.bed

# generate matched random control datasets

# identify all RE instances present in the reference genome
echo -e ">AvrII\nCCTAGG" > AvrII.fa
oligoMatch AvrII.fa ~/references/human/hg19.fa stdout \
| awk -v OFS="\t" '{split($4,enz,"+"); print $1,$2,$3,enz[1],$5,$6}' \
| sort -k1,1 -k2,2n \
> AvrII.hg19.bed

echo -e ">SpeI\nACTAGT" > SpeI.fa
oligoMatch SpeI.fa ~/references/human/hg19.fa stdout \
| awk -v OFS="\t" '{split($4,enz,"+"); print $1,$2,$3,enz[1],$5,$6}' \
| sort -k1,1 -k2,2n \
> SpeI.hg19.bed

echo -e ">NheI\nGCTAGC" > NheI.fa
oligoMatch NheI.fa ~/references/human/hg19.fa stdout \
| awk -v OFS="\t" '{split($4,enz,"+"); print $1,$2,$3,enz[1],$5,$6}' \
| sort -k1,1 -k2,2n \
> NheI.hg19.bed

echo -e ">MseI\nTTAA" > MseI.fa
oligoMatch MseI.fa ~/references/human/hg19.fa stdout \
| awk -v OFS="\t" '{split($4,enz,"+"); print $1,$2,$3,enz[1],$5,$6}' \
| sort -k1,1 -k2,2n \
> MseI.hg19.bed

# generate a bed file with all AvrII, NheI and SpeI site (digestion by enzyme mix to prepare AvrII libraries)
cat AvrII.hg19.bed NheI.hg19.bed SpeI.hg19.bed \
| sort -k1,1 -k2,2n \
> avrii_mix.hg19.bed

# generate multiple (here nb_mrc_per_is=10) matched random controls (mrc) for each insertion site (is)
# here restriction enzyme (re) used to generate the libraries were a pool of 3 re (AvrII, NheI, SpeI)

# calculate the distance of each is to the closest corresponding re site (non overlapping)
bedtools closest -D a -io -a ~/projects/iss/is_datasets/hiv/hiv.is_jurkat.avrii.hg19.bed -b avrii_mix.hg19.bed \
| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$(NF-3),$NF}' \
| sort -k1,1 -k2,2n \
> avrii_mix_is_dist.hg19.bed

bedtools closest -D ref -t first -iu -io -a avrii_mix.hg19.bed -b avrii_mix.hg19.bed \
> avrii_mix_downstream.tmp

bedtools closest -D ref -t first -id -io -a avrii_mix.hg19.bed -b avrii_mix.hg19.bed \
> avrii_mix_upstream.tmp

paste -d "\t" avrii_mix_upstream.tmp avrii_mix_downstream.tmp \
| awk -v OFS="\t" '$7!="." && $20!="." {print $1,$2,$3,$4,$5,$6,$13,$26}' \
> avrii_mix_dist.bed

# this step is long for AvrII mix (7h on professorx)
bash 170524_mrc_generator_simple_avrii.sh

# test
# head -20 avrii_mix_tmp.bed | awk -v OFS="\t" '{print $0,NR}' | sort -k1,1 -k2,2n > test_mrc
# bedtools closest -D ref -io -t first -a test_mrc -b avrii_mix.hg19.bed | sort -k7,7n | head -20
# head -20 avrii_mix_is_dist.hg19.bed


# generate multiple (here nb_mrc_per_is=10) matched random controls (mrc) for each insertion site (is)
# here restriction enzyme (re) used to generate the libraries was MseI

# calculate the distance of each is to the closest corresponding re site (non overlapping)
bedtools closest -D a -t first -io -a ~/projects/iss/is_datasets/hiv/hiv.is_jurkat.msei.hg19.bed -b MseI.hg19.bed \
| awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$(NF-3),$NF}' \
| sort -k1,1 -k2,2n \
> msei_is_dist.hg19.bed

bedtools closest -D ref -t first -iu -io -a MseI.hg19.bed -b MseI.hg19.bed \
> msei_downstream.tmp

bedtools closest -D ref -t first -id -io -a MseI.hg19.bed -b MseI.hg19.bed \
> msei_upstream.tmp

paste -d "\t" msei_upstream.tmp msei_downstream.tmp \
| awk -v OFS="\t" '$7!="." && $20!="." {print $1,$2,$3,$4,$5,$6,$13,$26}' \
> msei_dist.bed

# this step is very very long for MseI (
bash 170524_mrc_generator_simple_msei.sh

# test
# head -20 msei_tmp.bed | awk -v OFS="\t" '{print $0,NR}' | sort -k1,1 -k2,2n > test_mrc
# bedtools closest -D ref -io -t first -a test_mrc -b MseI.hg19.bed | sort -k7,7n | head -20
# head -20 msei_is_dist.hg19.bed

