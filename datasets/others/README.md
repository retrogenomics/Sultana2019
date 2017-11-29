# Origin and processing of published datasets

## Common processing
All datasets were:
- reformated as 2 bp-interval `.bed` files
- intersected with [hg19.helaAllowedGenomeSpace.bed](https://github.com/retrogenomics/iss/tree/master/annotations). In other words, it excludes DAC blacklisted regions from ENCODE and the Duke excluded regions, as well as unassembled contigs, gaps, alternate assemblies, and the Y chromosome.

## Naming conventions
We adopted a common naming nomenclature for files and individual insertion names.

### File names
Rule: `<genome version>.<mobile element>.<fragmentation method>.<insertion or locus>.<cell type>.bed`
Example: `hg19.hiv.msei.ins.jurkat.bed`

### Insertion names
Rule: `mobile element|cell type|fragmentation method|ins or loc|unique id`
Example: `hiv|jurkat|avrii|ins|0001`

## Sleeping Beauty and PiggyBac
Sleeping beauty (sb) and PiggyBac datasets were published in [Gogol-Doring A et al. Mol Ther. 2016](https://www.ncbi.nlm.nih.gov/pubmed/26755332). The datasets (GSE58744) were downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE58744&format=file). We only kept libraries prepared from sonicated DNA. The experiments were unstranded, but we assigned arbitrarily a positive orientation to comply with BED6 format specifications.
