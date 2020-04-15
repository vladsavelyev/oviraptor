#!/usr/bin/env python
from Bio import SeqIO
import sys
import os
from numpy import mean
from os.path import basename, dirname, join

fq1 = sys.argv[1]
fq2 = fq1.replace('R1', 'R2')
filtered_fq1 = join(dirname(fq1), 'filt_short', basename(fq1))
filtered_fq2 = join(dirname(fq1), 'filt_short', basename(fq2))

fq1_i = SeqIO.parse(open(fq1), "fastq")
fq2_i = SeqIO.parse(open(fq2), "fastq")

def _failed(rec):
    quals = rec.letter_annotations["phred_quality"]
    return len(rec.seq) < 125 or any(q < 10 for q in quals) or mean(quals) < 25

i = 0
written_i = 0
with open(filtered_fq1, 'w') as fq1_o, open(filtered_fq2, 'w') as fq2_o:
    for rec1, rec2 in zip(fq1_i, fq2_i):
        i += 1
        if not _failed(rec1) and not _failed(rec2):
            SeqIO.write(rec1, fq1_o, 'fastq')
            SeqIO.write(rec2, fq2_o, 'fastq')
            written_i += 1
        if i % 100000 == 0:
            print(f'Processed {i}, written {written_i}')

print('Done.')
print(f'Written {written_i}, filtered {i - written_i}')
