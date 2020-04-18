# Detection of oncoviruses in tumor whole genome samples

## Installation

Via conda:

```
conda create -n oncoviruses \
    -c vladsaveliev -c conda-forge -c bioconda -c defaults \
    python==3.7.3 bwa samtools sambamba snakemake ngs_utils click
conda activate oncoviruses
git clone --recursive git@github.com:umccr/oncoviruses.git
pip install -e oncoviruses
```

Installing conda is optional if you have the following tools installed and available in $PATH: python3, bwa, samtools, sambamba, snakemake

Usage:

```
oncoviruses tumor.bam -o results --genomes [path to umccrise genomes dir]
```

To use 10 CPUs:

```
oncoviruses tumor.bam -o results -t10
```

To analyse specifically HPV18:

```
oncoviruses tumor.bam -o results -v HPV18
```


## Preparing reference data

Creating combined reference:

```
samtools faidx genomes/hg38/hg38.fa \
    chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM \
    > genomes/hg38/hg38noalt.fa
cat genomes/viral/gdc-viral.fa >> genomes/hg38/hg38noalt_plus_gdcviral.fa
samtools faidx genomes/hg38/hg38noalt_plus_gdcviral.fa
bwa index genomes/hg38/hg38noalt_plus_gdcviral.fa
```












