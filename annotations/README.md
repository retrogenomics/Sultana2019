######### Downloaded

# hg19.wgEncodeDacMapabilityConsensusExcludable.bed
# The DAC Blacklisted Regions aim to identify a comprehensive set of regions in the human genome that have anomalous, unstructured, high signal/read counts in next gen sequencing experiments independent of cell line and type of experiment. There were 80 open chromatin tracks (DNase and FAIRE datasets) and 20 ChIP-seq input/control tracks spanning ~60 human tissue types/cell lines in total used to identify these regions with signal artifacts. These regions tend to have a very high ratio of multi-mapping to unique mapping reads and high variance in mappability. Some of these regions overlap pathological repeat elements such as satellite, centromeric and telomeric repeats. However, simple mappability based filters do not account for most of these regions. Hence, it is recommended to use this blacklist alongside mappability filters. The DAC Blacklisted Regions track was generated for the ENCODE project.

Downloaded as .bed file from UCSC Table Browser:
	[Assembly] hg19 >
	[Group] Mapping and Sequencing >
	[Track] Mappability >
	[Table] DAC Blacklist

# hg19.wgEncodeDukeMapabilityRegionsExcludable.bed
# The Duke Excluded Regions track displays genomic regions for which mapped sequence tags were filtered out before signal generation and peak calling for Open Chromatin: DNaseI HS and FAIRE tracks. This track contains problematic regions for short sequence tag signal detection (such as satellites and rRNA genes). The Duke Excluded Regions track was generated for the ENCODE project.

Downloaded as .bed file from UCSC Table Browser:
	[Assembly] hg19 >
	[Group] Mapping and Sequencing >
	[Track] Mappability >
	[Table] Duke Excluded


# hg19.gapsUnsorted.bed

Downloaded as .bed file from UCSC Table Browser:
	[Assembly] hg19 >
	[Group] Mapping and Sequencing >
	[Track] Gap >
	[Table] gap



######### Home-generated

# hg19.helaMainChromosomes.genome
> awk '$1~/^chr[0-9]+$|^chrX$/' hg19.genome \
> | sort -k1,1 \
> > hg19.helaMainChromosomes.genome

# hg19.helaMainChromosomes.bed
awk '$1~/^chr[0-9]+$|^chrX$/' hg19.bed \
| sort -k1,1 -k2,2n \
> hg19.helaMainChromosomes.bed

# hg19.gapsSorted.bed
sort -k1,1 -k2,2n hg19.gapsUnsorted.bed \
> hg19.gapsSorted.bed

# hg19.helaGapLessGenome.bed
bedtools subtract -a hg19.helaMainChromosomes.bed -b hg19.gapsSorted.bed \
| sort -k1,1 -k2,2n \
> hg19.helaGapLessGenome.bed
