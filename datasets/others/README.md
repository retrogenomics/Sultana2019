# Origin and processing of published datasets

## Common processing steps
All datasets were:
- reformated as 2 bp-interval `.bed` files
- intersected with [hg19.helaAllowedGenomeSpace.bed](https://github.com/retrogenomics/iss/tree/master/annotations). In other words, we excluded DAC blacklisted and Duke excluded regions, as well as unassembled contigs, gaps, alternate assemblies, and the Y chromosome.

## Naming conventions
We adopted a common naming nomenclature for files and individual insertion names.

### File names
- Rule: `<genome version>.<mobile element>.<fragmentation method>.<insertion or locus>.<cell type>.bed`
- Example: `hg19.hiv.msei.ins.jurkat.bed`

### Insertion names
- Rule: `mobile element|cell type|fragmentation method|ins or loc|unique id`
- Example: `hiv|jurkat|avrii|ins|0001`

## Sleeping Beauty and PiggyBac
Sleeping beauty (sb) and PiggyBac datasets were published in [Gogol-Doring A et al. Mol Ther. 2016](https://www.ncbi.nlm.nih.gov/pubmed/26755332). The datasets were downloaded from GEO ([GSE58744](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE58744&format=file)). We only kept libraries prepared from sonicated DNA. The experiments were unstranded, but we arbitrarily assigned a positive orientation to comply with BED6 format specifications. Genomic coordinates were converted by LiftOver from hg18 to hg19.

## Murine Leukemia Virus
Murine Leukemia Virus datasets were published in [LaFave MC et al. Nucleic Acids Res. 2014](https://www.ncbi.nlm.nih.gov/pubmed/?term=24464997) and were downloaded from the NHRGI website, [here](https://research.nhgri.nih.gov/software/GeIST/download.shtml). Genomic coordinates were already in hg19.

## Human Immunodefiency Virus
Human Immunodefiency Virus datasets were published in [Wang GP et al. Genome Res. 2007](https://www.ncbi.nlm.nih.gov/pubmed/17545577) and were downloaded from the author website, [here](http://microb230.med.upenn.edu/ucsc/hiv.wig.bed) as a `.wig` file. The data were unstranded, but we arbitrarily assigned a positive orientation to comply with BED6 format specifications. Genomic coordinates were converted by LiftOver from hg17 to hg19.

## Endogenous L1 elements
Endogenous L1 elements present in HeLaS3 genome were published in [Philippe C et al. eLife. 2016](https://www.ncbi.nlm.nih.gov/pubmed/27016617) and were downloaded from [Supplementary file 3](https://doi.org/10.7554/eLife.13926.018).
Genomic coordinates were already in hg19.
