# iss_scripts
Scripts used to obtain data for Sultana <I>et al.</I> manuscript.
- to call L1 insertions from sequencing data
- to format datasets
- to generate control datasets
- to generate figures

# Installation
## Dependencies
- [cutadapt](https://github.com/marcelm/cutadapt) (tested version: 1.14)
- [bwa](https://github.com/lh3/bwa) (tested version: 0.7.10-r789)
- [Picard tools](http://broadinstitute.github.io/picard/) (tested version: 1.136)
- [bedtools](https://github.com/arq5x/bedtools2) (tested version: 2.25.0)
- [seqtk](https://github.com/lh3/seqtk) (tested version: 1.0-r82-dirty; note that seqtk is only required if subsampling of sequencing data is used - to reduce time of analysis in tests)
- GNU grep/awk

## Procedure
