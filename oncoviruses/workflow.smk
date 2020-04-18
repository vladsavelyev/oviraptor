import sys
from os.path import isfile, join, basename, dirname, abspath
import csv
import yaml
from ngs_utils.file_utils import which, safe_mkdir
from ngs_utils.file_utils import get_ungz_gz
from ngs_utils.file_utils import splitext_plus
from ngs_utils.logger import warn
from hpc_utils import hpc
from oncoviruses import package_path


#############################
#### Reading parameters #####

INPUT_BAM = config['input_bam']
OUTPUT_DIR = abspath(config['output_dir'])
SAMPLE = config['sample_name']
VIRUS = config.get('virus', None)
if VIRUS:
    VIRUS = VIRUS.upper()
THREADS = config.get('cores', 10)

GENOME = config.get('genome', 'hg38')
if config.get('genomes_dir'):
    hpc.set_genomes_dir(config.get('genomes_dir'))

COMBINED_FA = config.get('combined_fa') or \
    hpc.get_ref_file(genome=GENOME, key='fa_plus_gdc_viruses', must_exist=False)
HOST_FA = config.get('host_fa') or \
    hpc.get_ref_file(genome=GENOME, key='fa', must_exist=False)
VIRUSES_FA = config.get('viruses_fa') or \
    hpc.get_ref_file(key='gdc_viral_fa', must_exist=False)
assert VIRUSES_FA and (COMBINED_FA or HOST_FA)

WORK_DIR = join(OUTPUT_DIR, 'work')


rule all:
    input: join(WORK_DIR, 'all.done')


rule extract_unmapped_and_mate_unmapped_reads:
    input:
        host_bam = INPUT_BAM
    output:
        host_bam_namesorted = join(WORK_DIR, f'step1_host_unmapped_or_mate_unmapped.namesorted.bam')
    threads: min(10, THREADS)
    shell:
         "sambamba view -t{threads} -fbam -F 'unmapped or mate_is_unmapped' {input.host_bam}"
         " | samtools sort -n -@{threads} -Obam -o {output.host_bam_namesorted}"


rule unmapped_and_mate_unmapped_reads_to_fastq:
    input:
        host_bam_namesorted = rules.extract_unmapped_and_mate_unmapped_reads.output.host_bam_namesorted,
    output:
        fq1 = join(WORK_DIR, f'step2_host_unmapped_or_mate_unmapped.R1.fq'),
        fq2 = join(WORK_DIR, f'step2_host_unmapped_or_mate_unmapped.R2.fq'),
    threads: max(10, THREADS)
    shell:
        "samtools fastq -@{threads} {input.host_bam_namesorted} -1 {output.fq1} -2 {output.fq2}"


if not VIRUS:
    # aligning to all viral genomes to figure out which one aligns the best
    rule bwa_unmapped_and_mate_unmapped_reads_to_gdc:
        input:
            fq1 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq1,
            fq2 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq2,
            gdc_fa = VIRUSES_FA,
        output:
            gdc_bam = join(WORK_DIR, 'detect_viral_reference', 'host_unmapped_or_mate_unmapped_to_gdc.bam')
        threads: THREADS
        shell:
            "bwa mem -Y -t{threads} {input.gdc_fa} {input.fq1} {input.fq2}"
            " | samtools sort -@{threads} -Obam -o {output.gdc_bam}"

    rule index_virus_bam:
        input:
            rules.bwa_unmapped_and_mate_unmapped_reads_to_gdc.output.gdc_bam,
        output:
            rules.bwa_unmapped_and_mate_unmapped_reads_to_gdc.output.gdc_bam + '.bai',
        shell:
            "samtools index {input}"

    mosdepth_work_dir = join(WORK_DIR, 'detect_viral_reference', 'mosdepth')
    mosdepth_prefix = join(mosdepth_work_dir, SAMPLE)
    rule mosdepth:
        input:
            bam = rules.bwa_unmapped_and_mate_unmapped_reads_to_gdc.output.gdc_bam,
            bai = rules.bwa_unmapped_and_mate_unmapped_reads_to_gdc.output.gdc_bam + '.bai',
            fai = VIRUSES_FA + '.fai',
        output:
            mosdepth_regions_bed_gz = mosdepth_prefix + '.regions.bed.gz',
            mosdepth_thresholds_bed_gz = mosdepth_prefix + '.thresholds.bed.gz',
        params:
            work_dir = mosdepth_work_dir,
            prefix = mosdepth_prefix,
        threads: THREADS
        run:
            awk_cmd = "awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}'"
            mosdepth_cmd = (
                f'mosdepth {params.prefix} {input.bam} -t{threads} -n --thresholds 1,5,25 '
                f'--by <({awk_cmd} {input.fai}) '
            )
            if hpc.name == 'vlad':
                # using dockerized mosdepth
                bam_dir = abspath(dirname(input.bam))
                fai_dir = abspath(dirname(input.fai))
                docker_mosdepth_cmd = (
                    f'docker run -v{bam_dir}:{bam_dir} -v{fai_dir}:{fai_dir} -v{params.work_dir}:/work '
                    f'quay.io/biocontainers/mosdepth:0.2.9--hbeb723e_0 ' +
                    mosdepth_cmd.replace(params.work_dir, '/work')
                )
                shell(docker_mosdepth_cmd)
            else:
                shell(mosdepth_cmd)

    ONCOVIRAL_SOURCE_URL = 'https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files'
    rule prioritize_viruses:
        input:
            mosdepth_regions_bed_gz = rules.mosdepth.output.mosdepth_regions_bed_gz,
            mosdepth_thresholds_bed_gz = rules.mosdepth.output.mosdepth_thresholds_bed_gz,
        output:
            oncoviruses_tsv = join(OUTPUT_DIR, 'prioritized_oncoviruses.tsv'),
        shell:
            "echo '## Viral sequences (from {ONCOVIRAL_SOURCE_URL}) found in unmapped reads' > {output.oncoviruses_tsv} &&"
            "echo '## Sample: {SAMPLE}' >> {output.oncoviruses_tsv} && "
            "echo '#virus\tsize\tdepth\t1x\t5x\t25x' >> {output.oncoviruses_tsv} && "
            "paste <(zcat {input.mosdepth_regions_bed_gz}) <(zgrep -v ^# {input.mosdepth_thresholds_bed_gz}) | "
            "awk 'BEGIN {{FS=\"\\t\"}} {{ print $1 FS $3 FS $4 FS $10/$3 FS $11/$3 FS $12/$3}}' | "
            "sort -n -r -k 5,5 >> {output.oncoviruses_tsv}"

    checkpoint select_viruses:
        input:
            tsv = rules.prioritize_viruses.output.oncoviruses_tsv,
            gdc_fa = VIRUSES_FA,
        output:
            selected_viruses_tsv = join(WORK_DIR, 'detect_viral_reference', 'integrated_viruses.tsv'),
            # virus_fas_dir = directory(join(WORK_DIR, 'viral_fa')),
        run:
            viruses = []
            with open(input.tsv) as f:
                for l in f:
                    if not l.startswith('#'):
                        virus, size, depth, t1, t5, t25 = l.strip().split('\t')
                        # at least 50% of a virus covered with at least 5x reads:
                        if float(t5) > 0.5:
                            viruses.append(virus)
            if not viruses:
                warn(
                    'Not found any viruses with significant coverage. '
                    'You can explore the full list at {input.tsv} and rerun with '
                    'the --virus option explicitly.')
                sys.exit(0)
            else:
                print(f'Found viruses: {", ".join(viruses)}')
                with open(output.selected_viruses_tsv, 'w') as out:
                    for v in viruses:
                        out.write(v + '\n')
                    # shell(f"samtools faidx {input.gdc_fa} {v} > {join(output.virus_fas_dir, f'{v}.fa')}")

rule create_viral_reference:
    input:
        gdc_fa = VIRUSES_FA,
    output:
        virus_fa = join(WORK_DIR, '{virus}', '{virus}.fa'),
        virus_bwt = join(WORK_DIR, '{virus}', '{virus}.fa.bwt'),  # one of bwa index files
    # params:
    #     virus =
        # selected_viruses_tsv = join(WORK_DIR, 'integrated_viruses.tsv'),
    shell:
        # if isfile(params.selected_viruses_tsv):
        #     viruses = open(params.selected_viruses_tsv).readlines()
        "samtools faidx {input.gdc_fa} {wildcards.virus} > {output.virus_fa}"
        " && bwa index {output.virus_fa}"

# aligning to specific viral sequence
rule bwa_unmapped_and_mateunmapped_to_viral_ref:
    input:
        fq1 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq1,
        fq2 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq2,
        virus_bwa_prefix = rules.create_viral_reference.output.virus_fa,
    output:
        virus_bam_possorted = join(WORK_DIR, 'step3_host_unmapped_and_bridging_reads_to_{virus}.possorted.bam')
    threads: THREADS
    shell:
        # using the polyidus bwa command.
        # -T1 = minimum score to output [default 30]
        # -a  = output all alignments for SE or unpaired PE
        # -Y  = use soft clipping for supplementary alignments
        "bwa mem -a -Y -t{threads} {input.virus_bwa_prefix} {input.fq1} {input.fq2}"
        " | samtools sort -@{threads} -Obam -o {output.virus_bam_possorted}"


# Extracts reads with at least 50 softclipped base stretches.
# E.g. 50 bases are represented in CIGAR with 50S. So we should look for
# any 3 digit numbers, or 2 digit numbers starting with 5..9:
sambamba_softclip_expr = '(cigar =~ /\d{3,}S/ or cigar =~ /[56789]\dS/) and not supplementary'

# Removing reads that are not helpful: fully unmapped pairs. In other words, keeping
# mapped pairs and pairs with an unmapped mate that can bridge us to the host genome
rule extract_viral_and_bridging_reads:
    input:
        virus_bam_possorted = rules.bwa_unmapped_and_mateunmapped_to_viral_ref.output.virus_bam_possorted,
    output:
        virus_bam_possorted = join(WORK_DIR, 'step4_host_unmapped_and_bridging_reads_to_{virus}.only_bridging_reads.possorted.bam')
    threads: max(10, THREADS)
    shell:
         "sambamba view -t{threads} -fbam -F 'not unmapped or not mate_is_unmapped or {sambamba_softclip_expr}'"
         " {input.virus_bam_possorted}"
         " | samtools sort -@{threads} -Obam -o {output.virus_bam_possorted}"

rule namesort_viral_bridging_bam:
    input:
        virus_bam_possorted = rules.extract_viral_and_bridging_reads.output.virus_bam_possorted,
    output:
        virus_bam_namesorted = join(WORK_DIR, 'step5_host_unmapped_and_bridging_reads_to_{virus}.only_bridging_reads.namesorted.bam')
    threads: THREADS
    shell:
        "samtools sort -n -@{threads} {input.virus_bam_possorted} -Obam -o {output.virus_bam_namesorted}"

# now aligning back to host genome
rule viral_bridging_reads_to_fastq:
    input:
        virus_bam_namesorted = rules.namesort_viral_bridging_bam.output.virus_bam_namesorted,
    output:
        fq1 = join(WORK_DIR, 'step6_{virus}_bridging.R1.fq'),
        fq2 = join(WORK_DIR, 'step6_{virus}_bridging.R2.fq'),
    threads: max(10, THREADS)
    shell:
        "samtools fastq -@{threads} {input.virus_bam_namesorted} -1 {output.fq1} -2 {output.fq2}"

rule create_combined_reference:
    input:
        host_fa = HOST_FA,
        viruses_fa = VIRUSES_FA,
    output:
        combined_fa = join(OUTPUT_DIR, 'combined_reference', 'host_plus_viruses.fa'),
        combined_bwt = join(OUTPUT_DIR, 'combined_reference', 'host_plus_viruses.fa.bwt'),
    shell:
        "cat {input.host_fa} {input.viruses_fa} > {output.combined_fa} && "
        "samtools faidx {output.combined_fa} && "
        "bwa index {output.combined_fa}"

# aligning to specific viral sequence
rule bwa_viral_bridging_to_comb_ref:
    input:
        fq1 = rules.viral_bridging_reads_to_fastq.output.fq1,
        fq2 = rules.viral_bridging_reads_to_fastq.output.fq2,
        bwa_prefix = COMBINED_FA if COMBINED_FA is not None else rules.create_combined_reference.output.combined_fa,
    output:
        comb_bam_possorted = join(WORK_DIR, 'step7_{virus}_bridging_to_comb_ref.possorted.bam')
    threads: THREADS
    shell:
        # using the polyidus bwa command.
        # -T1 = minimum score to output [default 30]
        # -a  = output all alignments for SE or unpaired PE
        # -Y  = use soft clipping for supplementary alignments
        "bwa mem -a -Y -t{threads} {input.bwa_prefix} {input.fq1} {input.fq2}"
        " | samtools sort -@{threads} -Obam -o {output.comb_bam_possorted}"

# # Removing reads that are not helpful: fully unmapped pairs. In other words, keeping
# # mapped pairs and pairs with an unmapped mate that can bridge us to the host genome
# rule extract_viral_and_bridging_reads:
#     input:
#         comb_bam_possorted = rules.bwa_viral_bridging_reads_to_comb_ref.output.comb_bam_possorted,
#     output:
#         comb_bam_possorted = join(WORK_DIR, 'step4_comb_ref_mapped_and_bridging.possorted.bam')
#     threads: THREADS
#     shell:
#          "sambamba view -t{threads} -fbam -F 'not unmapped or not mate_is_unmapped' {input.comb_bam_possorted}"
#          " | samtools sort -@{threads} -Obam -o {output.comb_bam_possorted}"

rule index_comb_ref_bam:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted
    output:
        bai = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted + '.bai'
    shell:
        "samtools index {input.bam}"

rule prep_manta:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
        bai = rules.index_comb_ref_bam.output.bai,
        ref = COMBINED_FA if COMBINED_FA is not None else rules.create_combined_reference.output.combined_fa,
    output:
        run_script = join(WORK_DIR, 'step8_{virus}_manta/runWorkflow.py'),
    params:
        manta_dir = join(WORK_DIR, 'step8_{virus}_manta')
    group: 'manta'
    shell:
        "configManta.py "
        "--bam={input.bam} "
        "--referenceFasta={input.ref} "
        "--exome "
        "--runDir={params.manta_dir}"

rule run_manta:
    input:
        run_script = rules.prep_manta.output.run_script,
    output:
        vcf = join(WORK_DIR, 'step8_{virus}_manta/results/variants/candidateSV.vcf.gz'),
    threads: max(10, THREADS)
    group: 'manta'
    shell:
        "{input.run_script} -m local"

rule filter_manta:
    input:
        vcf = rules.run_manta.output.vcf
    output:
        vcf = join(OUTPUT_DIR, 'breakpoints_{virus}.vcf.gz')
    shell:
        "bcftools view {input.vcf} | "
        "grep -P '^#|{wildcards.virus}' | "
        "bcftools filter -i \"IMPRECISE=0\" -Oz -o {output.vcf}"


# rule bwa_bridging_reads_to_genome:
#     input:
#         fq1 = rules.viral_bridging_reads_to_fastq.output.fq1,
#         fq2 = rules.viral_bridging_reads_to_fastq.output.fq2,
#     output:
#         human_bam_namesorted = join(WORK_DIR, '{virus}_bridging_reads_to_host.namesorted.bam')
#     params:
#         human_bwa_prefix = join(hpc.get_ref_file(genome=GENOME, key='bwa'), f'{GENOME}.fa')
#     threads: THREADS
#     shell:
#         # using the polyidus bwa command.
#         # -T1 = minimum score to output [default 30]
#         # -a  = output all alignments for SE or unpaired PE
#         # -C  = append FASTA/FASTQ comment to SAM output
#         # -Y  = use soft clipping for supplementary alignments
#         "bwa mem -T10 -t{threads} -a -C -Y {params.human_bwa_prefix} {input.fq1} {input.fq2}"
#         " | samtools sort -n -@{threads} -Obam -o {output.human_bam_namesorted}"
#
# rule bwa_bridging_reads_to_genome_index:
#     input:
#          human_bam_namesorted = rules.bwa_bridging_reads_to_genome.output.human_bam_namesorted,
#     output:
#          human_bam_possorted = join(WORK_DIR, '{virus}_bridging_reads_to_host.possorted.bam'),
#          human_bai_possorted = join(WORK_DIR, '{virus}_bridging_reads_to_host.possorted.bam.bai'),
#     threads: THREADS
#     shell:
#          'samtools sort {input.human_bam_namesorted} -@{threads} -Obam -o {output.human_bam_possorted}'
#          ' && samtools index {output.human_bam_possorted}'
#
# rule find_approximate_integrations:
#     input:
#         human_bam_namesorted = rules.bwa_bridging_reads_to_genome.output.human_bam_namesorted,
#         virus_bam_namesorted = rules.namesort_viral_bridging_bam.output.virus_bam_namesorted,
#     output:
#         integrations_tsv = join(OUTPUT_DIR, '{virus}_approx_integration_info.tsv')
#     shell:
#         'polyidus_find_approx_integrations.py {input.human_bam_namesorted} {input.virus_bam_namesorted} '
#         '-v {wildcards.virus} > {output.integrations_tsv}'
#
# rule find_exact_integrations:
#     input:
#         human_bam_possorted = rules.bwa_bridging_reads_to_genome_index.output.human_bam_possorted,
#         virus_bam_namesorted = rules.namesort_viral_bridging_bam.output.virus_bam_namesorted,
#         approx_integrations = rules.find_approximate_integrations.output.integrations_tsv,
#     output:
#         integrations_tsv = join(OUTPUT_DIR, '{virus}_exact_integrations.tsv')
#     shell:
#         'polyidus_find_exact_integrations.py {input.human_bam_possorted} {input.virus_bam_namesorted} '
#         '{input.approx_integrations} -v {wildcards.virus} > {output.integrations_tsv}'
#
# def aggregate_viruses_input_fn(wildcards):
#     if not VIRUS:
#         # # calling just to wait for the execution of the checkpoint rule:
#         # virus_fas_dir = checkpoints.select_viruses.get(**wildcards).output.virus_fas_dir
#         # # this function will check files existing with this name pattern, and get all wildcards matching these files
#         # virus = glob_wildcards(join(WORK_DIR, 'viral_fa', '{virus}.fa')).virus
#         selected_viruses_tsv = checkpoints.select_viruses.get(**wildcards).output.selected_viruses_tsv
#         viruses = [v.strip() for v in open(selected_viruses_tsv).readlines() if v.strip()]
#     else:
#         viruses = [VIRUS]
#
#     return expand(rules.find_exact_integrations.output.integrations_tsv, virus=viruses)

def aggregate_viruses_input_fn(wildcards):
    if not VIRUS:
        selected_viruses_tsv = checkpoints.select_viruses.get(**wildcards).output.selected_viruses_tsv
        viruses = [v.strip() for v in open(selected_viruses_tsv).readlines() if v.strip()]
    else:
        viruses = [VIRUS]
    return expand(rules.filter_manta.output.vcf, virus=viruses)

rule aggregate_over_viruses:
    input:
        aggregate_viruses_input_fn
        # rules.filter_manta.output.vcf
    output:
        join(WORK_DIR, 'all.done')
    shell:
        'touch {output}'



