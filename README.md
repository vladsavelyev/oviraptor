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














