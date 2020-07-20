# Detection of oncoviruses in tumor whole genome samples

```
oncoviruses tumor.bam -o results --host-fa hg38.fa
```

## Installation

Via conda:

```
conda create -n oncoviruses \
    -c vladsaveliev -c conda-forge -c bioconda -c defaults \
    python==3.7.3 minimap2 samtools sambamba bcftools mosdepth snakemake ngs_utils click
conda activate oncoviruses
git clone --recursive git@github.com:umccr/oncoviruses.git
pip install -e oncoviruses
```

Installing conda is optional if you have the following tools installed and available in $PATH:

- `python3`
- `minimap2`
- `samtools`
- `sambamba`
- `bcftools`
- `mosdepth`
- `snakemake`
- `lumpy` (on macOS, you can pull a dockerized version 
           with `docker pull quay.io/biocontainers/lumpy-sv:0.3.0--h0b85cd1_0`)

## Usage

The tool requires a whole genome seqeuncing BAM file as an input data, and the host (human) reference genome fasta file (ideally, hg38), which can be provided with `--host-fa` as follows:

```
oncoviruses tumor.bam -o results --host-fa hg38.fa
```

The tool is using an embedded ensemble annotation file to annotate the breakpoints with genes. You can override it with your own annotation file with `--host-gtf`, e.g.:

```
oncoviruses tumor.bam -o results --host-fa hg38.fa --host-gtf Homo_sapiens.GRCh38.gtf.gz
```

If the `--host-fa` is not provided, the tool will attempt to download into the output folder, which might take a while and around 3G of space.

The tool can make use of multiple cores. To use 10 CPUs:

```
oncoviruses tumor.bam -o results -t10
```

If you already know your candidate virus and just want to find the integration sites, use:

```
oncoviruses tumor.bam -o results -v HPV18
```

If you don't need integration sites and just want to find viral content, use:

```
oncoviruses tumor.bam -o results --only-detect
```










