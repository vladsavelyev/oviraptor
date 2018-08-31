#!/usr/bin/env python3
import sys
import os
from numpy import mean
import pysam


""" Usage:
./filter_unmapped_bam.py input.bam output.bam --10x
"""


bam = sys.argv[1]
bamf_bx_hqual_out = sys.argv[2]

# print('Loading whitelisted barcodes')
# with open('/g/data3/gx8/projects/Saveliev_10X/NA12878-10x/insertions/unmapped_or_mate_is_unmapped/4M-with-alts-february-2016.txt') as f:
#     barcodes = set(l.strip() for l in f.readlines())
# print(f'Read {len(barcodes)} barcodes: ')
# print('...')
# print()
# DON'T NEED TO CHECK BACODES BECAUSE THE WHILE LIST WAS CHECKED ALREADY IN PRE-PROCESSING

bamf = pysam.AlignmentFile(bam, "rb")
bamf_bx_hqual_f = pysam.AlignmentFile(bamf_bx_hqual_out, "wb", template=bamf)


writing_mates = []


def is_good_mapped_read(r_):
    """ We filter reads unmapped to a human genome, or having a mapped mate. 
        We want that mate to be a primary alignment, high quality and not duplicate. """
    return \
        not mate.is_unmapped and \
        not r_.is_secondary and \
        not r_.is_duplicate and \
        r_.mapping_quality != 0 and \
        r_.query_alignment_length >= 60


def is_good_read(r_):
    if not r_.is_paired: return False

    """ We also want the reads to have a high sequencing quality"""
    if r_.query_length < 125: return False
    if r_.is_qcfail: return False
    if not r_.query_qualities: return False
    hqual = all(q >= 10 for q in r_.query_qualities) and mean(r_.query_qualities) >= 25
    if not hqual: return False

    if len(sys.argv) > 3 and sys.argv[3] == '--10x':
        bx = r_.has_tag('BX')
        if not bx: return False

    return True


def resolve_mates(mates):
    # We want all high phred quaality pairs, where at least one read is mapped with a good quality:
    if any(is_good_mapped_read(mate) for mate in mates) and \
       all(is_good_read(mate) for mate in mates):
        for mate in mates:
            bamf_bx_hqual_f.write(mate)


for read in bamf:
    if not writing_mates or writing_mates[0].query_name == read.query_name:
        writing_mates.append(read)
    else:
        resolve_mates(writing_mates)
        writing_mates = [read]

    # TODO: try to reconstruct the original BX?


bamf.close()
bamf_bx_hqual_f.close()

# with open(filtered_fq1, 'w') as fq1_o, open(filtered_fq2, 'w') as fq2_o:
#     for rec1, rec2 in zip(fq1_i, fq2_i):
#         i += 1
#         if not _failed(rec1) and not _failed(rec2):
#             SeqIO.write(rec1, fq1_o, 'fastq')
#             SeqIO.write(rec2, fq2_o, 'fastq')
#             written_i += 1
#         if i % 100000 == 0:
#             print(f'Processed {i}, written {written_i}')

print('Done.')

