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
4. Process sequencing data and get L1 insertions by running `atlas-neo_all_samples.sh `:
```bash
cd iss/scripts
./atlas-neo_all_samples.sh 
```
5. Reformat L1 data and generate control files by running `prepare_L1neo_datasets.sh`
```bash
cd iss/datasets/l1neo
./prepare_L1neo_datasets.sh
```

