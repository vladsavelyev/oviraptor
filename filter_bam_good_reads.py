#!/usr/bin/env python3
import sys
import os
from numpy import mean
from os.path import basename, dirname, join
import pysam

""" Extract all good reads:
    - at least 125bp
    - with all bases quality >=10
    - with mean bases quality >=25
    - paired

Usage:
./filter_good_reads.py input.bam
"""

bam = sys.argv[1]
bamf_hqual_out = join(dirname(bam), 'lng_hqual', basename(bam))

bamf = pysam.AlignmentFile(bam, "rb")
bamf_hqual_f = pysam.AlignmentFile(bamf_hqual_out, "wb", template=bamf)

total_reads = 0
total_good = 0

for read in bamf:
    total_reads += 1

    lng = read.query_length >= 125
    if not read.query_qualities:
        continue

    hqual = all(q >= 10 for q in read.query_qualities) and mean(read.query_qualities) >= 25

    if read.is_paired and lng and hqual:
        bamf_hqual_f.write(read)
        total_good += 1

bamf.close()
bamf_hqual_f.close()

print(f'Done. Written {total_good} out of {total_reads}')