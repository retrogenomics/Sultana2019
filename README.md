# iss
## Background
Scripts used to obtain data for Sultana <I>et al.</I> manuscript.
- to call L1 insertions from sequencing data
- to format datasets
- to generate control datasets
- to generate figures

## Installation
### Dependencies
- [cutadapt](https://github.com/marcelm/cutadapt) (tested version: 1.14)
- [bwa](https://github.com/lh3/bwa) (tested version: 0.7.16a)
- [Picard tools](http://broadinstitute.github.io/picard/) (tested version: 1.136)
- [bedtools](https://github.com/arq5x/bedtools2) (tested version: 2.25.0)
- [seqtk](https://github.com/lh3/seqtk) (tested version: 1.0; note that seqtk is only required if subsampling of sequencing data is used - to reduce time of analysis in tests)
- [GNU parallel](https://www.gnu.org/software/parallel/) (tested version: 20150522)
- GNU grep/awk

### Other requirements
- A human reference genome sequence (ex:`hg19.fa`), 
- and its bwa index (ex: `hg19.fa.amb, .ann, .bwt, .fai, .pac, .sa`)

### Procedure

1. Download from GitHub:
```
git clone https://github.com/retrogenomics/iss.git
```
2. Edit the `CONFIG` file in the `iss/scripts/` folder according to your configuration
3. Edit the experiment configuration file (sequencing metadata) according to your sequencing data.
   - This file is stored in the `iss/experiments/` folder and its name.
   - Its location should be included in the `CONFIG`file as `EXP_CONFIG="${EXPERIMENTS}/<your file name>"`.
   - `atlas-neo-R01-to-R09.txt` is the default example.

## How to use?
### To process sequencing data and get L1 insertions
```bash
cd iss/scripts
./atlas-neo_all_samples.sh
```

### To generate control L1 files

Go to the `(...)/iss/datasets/l1neo` directory and run the following commands for target loci (loc, non-redundant) and insertions (ins, can be redundant if several L1 insert at the same nucleotidic position).

Note that control datasets are generated using the same filters as the experimental dataset, i.e. using the `hg19.helaAllowedGenomeSpace.bed` file to keep only insertions which are in chr1-22 and chrX, and to exclude blacklisted regions (ENCODE, Duke).

1. Random datasets

```bash
../../scripts/random_generator.sh \
	-N 1000 \
	-n "l1neo|insilico|soni|random_loc" \
	-a ../../annotations/hg19.helaAllowedGenomeSpace.bed \
	-g ../../annotations/hg19.genome \
	-i hg19.l1neo.soni.loc.helas3.bed \
	-o random_loc \
	-p hg19.l1neo.soni.random_loc.
```

```bash	
../../scripts/random_generator.sh \
	-N 1000 \
	-n "l1neo|insilico|soni|random_ins" \
	-a ../../annotations/hg19.helaAllowedGenomeSpace.bed \
	-g ../../annotations/hg19.genome \
	-i hg19.l1neo.soni.ins.helas3.bed \
	-o random_ins \
	-p hg19.l1neo.soni.random_ins.
```

2. Matched random control (mrc) datasets

```bash
../../scripts/mrc_generator.sh \
	-N 1000 \
	-n "l1neo|insilico|soni|mrc_loc" \
	-a ../../annotations/hg19.helaAllowedGenomeSpace.bed \
	-g ../../annotations/hg19.genome \
	-f ~/references/human/hg19.fa \
	-i hg19.l1neo.soni.loc.helas3.bed \
	-o mrc_loc \
	-p hg19.l1neo.soni.mrc_loc.
```

```bash
../../scripts/mrc_generator.sh \
	-N 1000 \
	-n "l1neo|insilico|soni|mrc_ins" \
	-a ../../annotations/hg19.helaAllowedGenomeSpace.bed \
	-g ../../annotations/hg19.genome \
	-f ~/references/human/hg19.fa \
	-i hg19.l1neo.soni.ins.helas3.bed \
	-o mrc_ins \
	-p hg19.l1neo.soni.mrc_ins.
```
