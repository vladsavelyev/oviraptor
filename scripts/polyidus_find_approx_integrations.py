#!/usr/bin/env python
import click
from ngs_utils.utils import set_locale; set_locale()
import sys
from ngs_utils import logger
import polyidus.src.utils as polyidus


@click.command()
@click.argument('host_bam_namesorted', type=click.Path(exists=True))
@click.argument('virus_bam_namesorted', type=click.Path(exists=True))
@click.option('-v', '--virus', 'virus_name')
@click.option('-s', '--signle-end', 'signle_end', is_flag=True, default=False)
@click.option('-o', '--output-tsv', 'output_tsv', type=click.Path())

def main(host_bam_namesorted, virus_bam_namesorted, virus_name=None, signle_end=None,
         output_tsv=None):

    ad_list = polyidus.find_integration(host_bam_namesorted, virus_bam_namesorted,
                                        paired_end=not signle_end)
    if len(ad_list) == 0:
        logger.warn(f'Not found any candidate integration sites' +
                    (f' for virus {virus_name}' if virus_name else ''))
        sys.exit(0)

    header = [
        "HostFile", "ViralFile", "ChromHost",
        "PositionHost", "StrandHost",
        "ChromViral", "PositionViral", "StrandViral",
        "NumberReads", "Score",
        "Normalized.Score", "ReadNames"
    ]
    out = open(output_tsv, 'w') if output_tsv else sys.stdout
    for line in [header] + ad_list:
        out.write("\t".join(line) + "\n")


if __name__ == '__main__':
    main()
