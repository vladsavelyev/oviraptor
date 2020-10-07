#!/bin/bash
test -e chr8.fa.gz || wget --no-verbose https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr8.fa.gz
test -e chr8.fa.gz.gzi || gunzip chr8.fa.gz && bgzip chr8.fa && samtools faidx chr8.fa.gz
oviraptor tumor.bam -o results --host-fa chr8.fa.gz

