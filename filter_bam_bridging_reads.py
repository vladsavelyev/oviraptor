#!/usr/bin/env python3
import sys
import os
from numpy import mean
import pysam

""" Extract all read groups where all mates are good reads:
    - at least 125bp
    - with all bases quality >=10
    - with mean bases quality >=25
    - paired
    - not qc_failed
    And at least one mate is good mapped:
    - mapped
    - not secondary
    - not duplicated
    - mapping quality !=0
    - alignment length >=60

Usage:
# INPUT: name-sorted BAM file
./filter_good_discordant_mates.py input.name_sorted.bam output.bam --10x
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
        We want that mate to be a primary alignment, high quality and not duplicate.
    """
    return \
        not r_.is_unmapped and \
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
    """ We want all high phred quality pairs, that either all are unmapped, or at least one read is mapped with a good quality.
    """
    if not all(is_good_read(aln) for aln in mates):
        return []
    mapped_mates = [aln for aln in mates if not aln.is_unmapped]
    if not mapped_mates:
        return []
    if not all(is_good_mapped_read(aln) for aln in mapped_mates):
        return []
    return mates
    

for aln in bamf:
    if not writing_mates or writing_mates[0].query_name == aln.query_name:
        writing_mates.append(aln)
    else:
        for prev_aln in resolve_mates(writing_mates):
            bamf_bx_hqual_f.write(prev_aln)
        writing_mates = [aln]


bamf.close()
bamf_bx_hqual_f.close()