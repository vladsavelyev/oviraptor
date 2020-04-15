#!/usr/bin/env python
import click
from ngs_utils.utils import set_locale; set_locale()
import sys
from ngs_utils import logger
import numpy as np
import polyidus.src.utils as polyidus


@click.command()
@click.argument('host_bam_indexed', type=click.Path(exists=True))
@click.argument('viral_bam_namesorted', type=click.Path(exists=True))
@click.argument('approx_integrations_tsv', type=click.Path(exists=True))
@click.option('-v', '--virus', 'virus_name')
@click.option('-o', '--output-tsv', 'output_tsv', type=click.Path())

def main(host_bam_indexed, viral_bam_namesorted, approx_integrations_tsv,
         virus_name='virus', output_tsv=None):
    # Defaults
    minreads = 2
    minscore = 0.6
    chimer_df = polyidus.get_chimer(approx_integrations_tsv, minreads, minscore)

    ad_lists = []
    keys = ["ChromHost", "Starts", "Ends", "IntegrationHosts",
            "ChromVirus", "VirusStarts", "VirusEnds",
            "IntegrationVirus", "NumberReads"]
    for i in range(chimer_df.shape[0]):
        chromhost = str(chimer_df.iloc[i, 2])
        if "alt" not in chromhost and "rand" not in chromhost:
            readnames = chimer_df.iloc[i, -1].split(", ")
            start_host = int(
                np.mean([
                    int(each) for each in
                    chimer_df.iloc[i, 3].split(", ")
                ])
            )
            dictposes = polyidus.get_host_pos(
                viral_bam_namesorted, host_bam_indexed, chromhost,
                start_host, readnames)
            LEN_VIR = len(dictposes["ChromVirus"]) > 0
            if LEN_VIR and dictposes["ChromVirus"][0] != "NA":
                for j in range(len(dictposes["Starts"])):
                    ad_list = []
                    for key in keys:
                        ad_list.append(str(dictposes[key][j]))
                    ad_list.append(virus_name)
                    ad_lists.append(ad_list)

    if len(ad_lists) == 0:
        logger.warn(f'Could not detect any exact integration sites for the virus {virus_name}, '
             f'see approximate integrations at {approx_integrations_tsv}')
        sys.exit(0)

    header = ["Chrom", "Start", "End", "IntegrationSite",
              "ChromVirus", "VirusStart", "VirusEnd",
              "ViralIntegrationSite",
              "NumberOfSupportingReads", "SampleName"]

    out = open(output_tsv, 'w') if output_tsv else sys.stdout
    out.write('\t'.join(header) + '\n')
    for ad_list in ad_lists:
        out.write("\t".join(ad_list) + "\n")


if __name__ == '__main__':
    main()
