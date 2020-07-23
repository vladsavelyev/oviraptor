import os
from collections import defaultdict
from os.path import isfile, join, basename, dirname, abspath
import subprocess
from ngs_utils.call_process import run_simple
from ngs_utils.file_utils import open_gzipsafe, get_ungz_gz
from ngs_utils.logger import warn, critical
from ngs_utils.vcf_utils import count_vars
from oviraptor import package_path


#############################
#### Reading parameters #####
## Inputs
INPUT_BAM = config['input_bam']
OUTPUT_DIR = abspath(config['output_dir'])
RESULT_PATH = join(OUTPUT_DIR, 'breakpoints.vcf.gz')
PRIO_TSV = join(OUTPUT_DIR, 'prioritized_oncoviruses.tsv')

## Other parameters
SAMPLE = config['sample_name']
THREADS = config.get('cores', 10)
WORK_DIR = join(OUTPUT_DIR, 'work')

SV_CALLER = 'lumpy'
ALIGNER = config.get('aligner', 'minimap2')  # bwa

VIRUSES = config.get('viruses', None)
if VIRUSES:
    VIRUSES = [v.upper() for v in VIRUSES.split(',')]
ONLY_DETECT = config.get('only_detect', False)

## Reference files
GENOME = 'hg38'
COMBINED_FA = config.get('combined_fa')
HOST_FA = config.get('host_fa')
VIRUSES_FA = config.get('viruses_fa')
pyens = None
try:
    from reference_data import api as refdata
except:
    pass
else:
    if config.get('genomes_dir'):
        refdata.find_genomes_dir(config.get('genomes_dir'))
    COMBINED_FA = COMBINED_FA or refdata.get_ref_file(genome=GENOME, key='fa_plus_gdc_viruses', must_exist=False)
    HOST_FA     = HOST_FA     or refdata.get_ref_file(genome=GENOME, key='fa', must_exist=False)
    VIRUSES_FA  = VIRUSES_FA  or refdata.get_ref_file(genome=GENOME, key='gdc_viral_fa', must_exist=False)
VIRUSES_FA = VIRUSES_FA or join(package_path(), 'data', 'gdc-viral.fa')
assert HOST_FA or COMBINED_FA

# By default, using the built-in BED in package_path()/data/hg38_noalt.genes.bed.gz
# But if the user provides the GTF file, using it instead
GTF_PATH = config.get('gtf_file')

rule all:
    input:
        RESULT_PATH if not ONLY_DETECT else [],
        PRIO_TSV if ONLY_DETECT or not VIRUSES else [],


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
        single = join(WORK_DIR, f'step2_host_unmapped_or_mate_unmapped.single.fq'),
    threads: max(10, THREADS)
    shell:
        "samtools fastq -@{threads} {input.host_bam_namesorted} -1 {output.fq1} -2 {output.fq2} -s {output.single}"


if not VIRUSES:
    # aligning to all viral genomes to figure out which one aligns the best
    rule bwa_unmapped_and_mate_unmapped_reads_to_gdc:
        input:
            fq1 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq1,
            fq2 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq2,
            gdc_fa = VIRUSES_FA,
        output:
            gdc_bam = join(WORK_DIR, 'detect_viral_reference', 'host_unmapped_or_mate_unmapped_to_gdc.bam')
        threads: THREADS
        params:
            rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
        run:
            for i in range(0, 100):
                tmp_bam = output.gdc_bam + '.tmp.{:0>4d}.bam'.format(i)
                if isfile(tmp_bam):
                    os.system(f'rm {tmp_bam}')
                else:
                    break

            if ALIGNER == 'bwa':
                # -Y  = use soft clipping for supplementary alignments
                shell(f"bwa mem -Y -t{threads} -R '{params.rg}'"
                      f" {input.gdc_fa} {input.fq1} {input.fq2} "
                      f" | samtools sort -@{threads} -Obam -o {output.gdc_bam}")
            else:
                # -Y  = use soft clipping for supplementary alignments
                shell(f"minimap2 -ax sr -Y -t{threads} -R '{params.rg}'"
                      f" {input.gdc_fa} {input.fq1} {input.fq2}"
                      f" | samtools sort -@{threads} -Obam -o {output.gdc_bam}")

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
        shell:
            'mosdepth {params.prefix} {input.bam} -t{threads} -n --thresholds 1,5,25 '
            '--by <(awk \'BEGIN {{FS="\\t"}}; {{print $1 FS "0" FS $2}}\' {input.fai})'

    # we need at least one of these conditions to call significance:
    MIN_1x_PCT = 50.0  #  % of the viral sequence that must be covered at at least 1x (probably non-integrating)
    MIN_5x_LEN = 300   #  viral base pairs must be covered at at least 5x (which is amplified, thus integrating)
    ONCOVIRAL_SOURCE_URL = 'https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files'
    rule prioritize_viruses:
        input:
            md_regions_bed_gz = rules.mosdepth.output.mosdepth_regions_bed_gz,
            md_thresholds_bed_gz = rules.mosdepth.output.mosdepth_thresholds_bed_gz,
        output:
            tsv = PRIO_TSV,
        params:
            completeness_fraction =MIN_1x_PCT / 100.0
        shell:
            "echo '## Viral sequences (from {ONCOVIRAL_SOURCE_URL}) found in unmapped reads' > {output.tsv} &&"
            "echo '## Sample: {SAMPLE}' >> {output.tsv} && "
            "echo '## Minimal completeness: {MIN_1x_PCT}% at 1x or {MIN_5x_LEN}bp at 5x' >> {output.tsv} && "
            "echo '#virus\tsize\tdepth\t1x\t5x\t25x\tsignificance' >> {output.tsv} && "
            "paste <(gunzip -c {input.md_regions_bed_gz}) <(zgrep -v ^# {input.md_thresholds_bed_gz}) | "
            "awk 'BEGIN {{FS=\"\\t\"}} {{ printf(\"%s\\t%d\\t%3.1f\\t%3.3f\\t%3.3f\\t%3.3f\\t%s\\n\", "
            "$1, $3, $4, $9/$3, $10/$3, $11/$3, (($10>{MIN_5x_LEN} || $11/$3>{params.completeness_fraction}) "
            "? \"significant\" : \".\")) }}' | "
            "sort -n -r -k5,5 -k6,6 -k4,4 -k3,3 >> {output.tsv}"

    checkpoint select_viruses:
        input:
            tsv = rules.prioritize_viruses.output.tsv,
            gdc_fa = VIRUSES_FA,
        output:
            selected_viruses_tsv = join(WORK_DIR, 'detect_viral_reference', 'present_viruses.txt'),
        run:
            viruses = []
            with open(input.tsv) as f:
                for l in f:
                    if not l.startswith('#'):
                        virus, size, depth, t1, t5, t25, significance = l.strip().split('\t')
                        if significance != '.':
                            viruses.append(virus)
            if not viruses:
                warn(
                    f'Not found any viruses with significant coverage. '
                    f'You can explore the full list at {input.tsv} and rerun with '
                    f'the --virus option explicitly.')
                shell(f'touch {output.selected_viruses_tsv}')
            else:
                print(f'Found viruses: {", ".join(viruses)}')
                with open(output.selected_viruses_tsv, 'w') as out:
                    for v in viruses:
                        out.write(v + '\n')


rule create_viral_reference:
    input:
        gdc_fa = VIRUSES_FA,
    output:
        virus_fa = join(WORK_DIR, '{virus}', '{virus}.fa'),
    shell:
        "samtools faidx {input.gdc_fa} {wildcards.virus} > {output.virus_fa}"

rule index_viral_reference:
    input:
        join(WORK_DIR, '{virus}', '{virus}.fa'),
    output:
        join(WORK_DIR, '{virus}', '{virus}.fa.bwt'),  # one of bwa index files
    shell:
        "bwa index {input}"

# aligning to specific viral sequence
rule bwa_unmapped_and_mateunmapped_to_viral_ref:
    input:
        fq1 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq1,
        fq2 = rules.unmapped_and_mate_unmapped_reads_to_fastq.output.fq2,
        virus_fa = rules.create_viral_reference.output.virus_fa,
        virus_bwt = rules.index_viral_reference.output if ALIGNER == 'bwa' else [],
    output:
        virus_bam_possorted = join(WORK_DIR, 'step3_host_unmapped_and_bridging_reads_to_{virus}.possorted.bam')
    threads: THREADS
    params:
        rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
    run:
        if ALIGNER == 'bwa':
            # using the polyidus bwa command.
            # -T1 = minimum score to output [default 30]
            # -a  = output all alignments for SE or unpaired PE
            # -Y  = use soft clipping for supplementary alignments
            shell(f"bwa mem -a -Y -t{threads} -R '{params.rg}'"
                  f" {input.virus_fa} {input.fq1} {input.fq2}"
                  f" | samtools sort -@{threads} -Obam -o {output.virus_bam_possorted}")
        else:
            # -Y  = use soft clipping for supplementary alignments
            shell(f"minimap2 -ax sr -Y -t{threads} -R '{params.rg}'"
                  f" {input.virus_fa} {input.fq1} {input.fq2}"
                  " | samtools sort -@{threads} -Obam -o {output.virus_bam_possorted}")


# Extracts reads with at least 50 soft-clipped base stretches.
# E.g. 50 bases are represented in CIGAR with 50S. So we should look for
# any 3 digit numbers, or 2 digit numbers starting with 5..9:
sambamba_softclip_expr = '(cigar =~ /\d{3,}S/ or cigar =~ /[56789]\dS/) and not supplementary'

# Removing reads that are not helpful: fully unmapped pairs. In other words, keeping
# mapped pairs and pairs with an unmapped mate that can bridge us to the host genome
rule extract_viral_and_bridging_reads:
    input:
        virus_bam_possorted = rules.bwa_unmapped_and_mateunmapped_to_viral_ref.output.virus_bam_possorted,
    output:
        virus_bam_possorted = join(WORK_DIR,
           'step4_host_unmapped_and_bridging_reads_to_{virus}.only_bridging_reads.possorted.bam')
    threads: max(10, THREADS)
    shell:
         "sambamba view -t{threads} -fbam -F 'not unmapped or not mate_is_unmapped or {sambamba_softclip_expr}'"
         " {input.virus_bam_possorted}"
         " | samtools sort -@{threads} -Obam -o {output.virus_bam_possorted}"

rule namesort_viral_bridging_bam:
    input:
        virus_bam_possorted = rules.extract_viral_and_bridging_reads.output.virus_bam_possorted,
    output:
        virus_bam_namesorted = join(WORK_DIR,
            'step5_host_unmapped_and_bridging_reads_to_{virus}.only_bridging_reads.namesorted.bam')
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
        combined_fa   = COMBINED_FA or join(OUTPUT_DIR, 'combined_reference', 'host_plus_viruses.fa'),
        combined_fai = (COMBINED_FA or join(OUTPUT_DIR, 'combined_reference', 'host_plus_viruses.fa')) + '.fai',
    run:
        if input.host_fa.endswith('.gz'):
            shell(f"gunzip -c {input.host_fa} > {output.combined_fa}")
        else:
            shell(f"cat {input.host_fa} > {output.combined_fa}")
        if input.viruses_fa.endswith('.gz'):
            shell(f"gunzip -c {input.viruses_fa} >> {output.combined_fa}")
        else:
            shell(f"cat {input.viruses_fa} >> {output.combined_fa}")
        shell(f"samtools faidx {output.combined_fa}")

rule index_combined_reference:
    input:
        COMBINED_FA or join(OUTPUT_DIR, 'combined_reference', 'host_plus_viruses.fa'),
    output:
        (COMBINED_FA or join(OUTPUT_DIR, 'combined_reference', 'host_plus_viruses.fa')) + '.bwt',
    shell:
        "bwa index {input}"

# aligning to specific viral sequence
rule bwa_viral_bridging_to_comb_ref:
    input:
        fq1 = rules.viral_bridging_reads_to_fastq.output.fq1,
        fq2 = rules.viral_bridging_reads_to_fastq.output.fq2,
        host_fa  = rules.create_combined_reference.output.combined_fa,
        host_bwt = rules.index_combined_reference.output if ALIGNER == 'bwa' else [],
    output:
        comb_bam_possorted = join(WORK_DIR, 'step7_{virus}_bridging_to_comb_ref.possorted.bam')
    threads: THREADS
    params:
        rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
    run:
        if ALIGNER == 'bwa':
            # using the polyidus bwa command.
            # -T1 = minimum score to output [default 30]
            # -a  = output all alignments for SE or unpaired PE
            # -Y  = use soft clipping for supplementary alignments
            shell(f"bwa mem -a -Y -t{threads} -R '{params.rg}'"
                  f" {input.host_fa} {input.fq1} {input.fq2}"
                  f" | samtools sort -@{threads} -Obam -o {output.comb_bam_possorted}")
        else:
            # -Y  = use soft clipping for supplementary alignments
            shell(f"minimap2 -ax sr -Y -t{threads} -R '{params.rg}'"
                  f" {input.host_fa} {input.fq1} {input.fq2}"
                  f" | samtools sort -@{threads} -Obam -o {output.comb_bam_possorted}")

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

rule extract_discordant:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
    output:
        bam = join(WORK_DIR, 'step8_{virus}_lumpy/disc.bam')
    group: 'lumpy'
    shell:
        'samtools view -h -F 1294 {input.bam} | '
        'samtools sort -Obam -o {output.bam}'

rule extract_split:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
    output:
        bam = join(WORK_DIR, 'step8_{virus}_lumpy/split.bam')
    params:
        extractSplitReads = join(package_path(), 'lumpy', 'extractSplitReads_BwaMem')
    group: 'lumpy'
    shell:
        'samtools view -h {input.bam} '
        '| {params.extractSplitReads} -i stdin '
        '| samtools sort -Obam -o {output.bam}'

rule lumpy_histo_bam_subset:
    input:
        bam = INPUT_BAM,
    output:
        bam = join(WORK_DIR, 'step8_{virus}_lumpy/10kreads.bam')
    params:
        pairend_distro = join(package_path(), 'lumpy', 'pairend_distro.py')
    group: 'lumpy'
    run:
        import pysam
        inp_bam = pysam.AlignmentFile(input.bam, 'rb')
        out_bam = pysam.AlignmentFile(output.bam, 'wb', template=inp_bam)
        for i, read in enumerate(inp_bam):
            if i < 10000:
                out_bam.write(read)
            else:
                break

rule lumpy_histo:
    input:
        join(WORK_DIR, 'step8_{virus}_lumpy/10kreads.bam')
    output:
        join(WORK_DIR, 'step8_{virus}_lumpy/histo.txt')
    params:
        pairend_distro = join(package_path(), 'lumpy', 'pairend_distro.py')
    group: 'lumpy'
    shell:
        'samtools view {input} | '
        'python {params.pairend_distro} '
        '--read_length 151 '
        '-X 4 '
        '-N 10000 '
        '-o {output}'

rule run_lumpy:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
        disc = rules.extract_discordant.output.bam,
        split = rules.extract_split.output.bam,
        histo = rules.lumpy_histo.output,
        blacklist = join(package_path(), 'data', 'encode4_unified_blacklist.bed.gz'),
    output:
        vcf = join(WORK_DIR, 'step8_{virus}_lumpy.vcf'),
    params:
        lumpy = join(package_path(), 'lumpy', 'lumpy'),
        image = 'quay.io/biocontainers/lumpy-sv:0.3.0--h0b85cd1_0'
    group: 'lumpy'
    run:
        tool_cmd = (
            f'-mw 4 '
            f'-tt 0 '
            f'-x {input.blacklist} '
            f'-pe id:sample,bam_file:{input.disc},histo_file:{input.histo},mean:500,stdev:50,read_length:151,'
            f'min_non_overlap:151,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 '
            f'-sr id:sample,bam_file:{input.split},back_distance:10,weight:1,min_mapping_threshold:20 '
            f'> {output.vcf}'
        )
        try:
            run_simple(f'{params.lumpy} 2>&1 | grep -q Program')
        except:
            # otherwise, try one from docker (will pull automatically):
            try:
                run_simple(f'docker images -q {params.image} 2>/dev/null')
            except:
                critical(f'Can\'t run lumpy either from {params.lumpy}, or using Docker')
            else:
                volumes_dict = dict()
                for inp_path in input:
                    volumes_dict[dirname(inp_path)] = dirname(inp_path)
                volumes_arg = " ".join(f"-v{k}:{v}" for k, v in volumes_dict.items())
                shell(
                    f'docker run ' 
                    f'{volumes_arg} '
                    f'{params.image} '
                    f'bash -c "lumpy {tool_cmd}"'
                )
        else:
            shell(f'{params.lumpy} {tool_cmd}')

# - gsort is not working because it sorts the header and puts the first line
# in the VCF in the middle, leading to downstream failues
# 0 bcftools requires all contigs to be definded in the header and not working as well
rule sort_vcf:
    input:
        vcf = join(WORK_DIR, 'step8_{virus}_lumpy.vcf'),
        fai = rules.create_combined_reference.output.combined_fa + '.fai',
    output:
        vcf = join(WORK_DIR, 'step9_{virus}_lumpy.sorted.vcf'),
    run:
        sort_order = dict()
        with open(input.fai) as f:
            for i, l in enumerate(f):
                chrom = l.split('\t')[0]
                sort_order[chrom] = i
        hdr_lines = []
        rec_split_lines = []
        with open(input.vcf) as f:
            for l in f:
                if l.startswith('#'):
                    hdr_lines.append(l)
                else:
                    rec_split_lines.append(l.strip().split('\t'))
        with open(output.vcf, 'w') as f:
            for l in hdr_lines:
                f.write(l)
            for fields in sorted(rec_split_lines, key=lambda fs: (sort_order[fs[0]], int(fs[1]))):
                f.write('\t'.join(fields) + '\n')

rule run_manta:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
        bai = rules.index_comb_ref_bam.output.bai,
        ref = rules.create_combined_reference.output.combined_fa,
    output:
        vcf = join(WORK_DIR, 'step8_{virus}_manta/results/variants/candidateSV.vcf.gz'),
    params:
        work_dir = join(WORK_DIR, 'step8_{virus}_manta'),
        image = 'quay.io/biocontainers/manta:1.6.0--py27_0',
    group: 'manta'
    run:
        run_script = join(params.work_dir, "runWorkflow.py")
        shell(f'test -e {run_script} && rm {run_script} || exit 0')  # remove if exists
        tool_cmd = (
            f'configManta.py '
            f'--bam={input.bam} '
            f'--referenceFasta={input.ref} '
            f'--exome '
            f'--runDir={params.work_dir} && ls {params.work_dir} && '
            f'{run_script} -m local 2>&1'
            f' | grep -v "Adding command task"'
            f' | grep -v "Launching command task"'
            f' | grep -v "Task initiated on local node"'
            f' | grep -v "launched from master workflow"'
            f' 1>&2'
        )
        if subprocess.run(f'docker images -q {params.image} 2>/dev/null', shell=True).returncode == 0:
            bam_dir = abspath(dirname(input.bam))
            ref_dir = abspath(dirname(input.ref))
            shell(
                f'docker run ' 
                f'-v{bam_dir}:{bam_dir} ' 
                f'-v{ref_dir}:{ref_dir} ' 
                f'-v{params.work_dir}:/work '
                f'{params.image} '
                f'bash -c "{tool_cmd.replace(params.work_dir, "/work")}"'
            )
        else:
            PY2_ENV_PATH = config.get('py2_env_path')
            py2_conda_cmd = ''
            if PY2_ENV_PATH:
                py2_conda_cmd = 'export PATH=' + PY2_ENV_PATH + '/bin:$PATH; '
            shell(f'{py2_conda_cmd} {tool_cmd}')

if SV_CALLER == 'manta':
    sv_output_rule = rules.run_manta
else:
    sv_output_rule = rules.sort_vcf

rule filter_vcf:
    input:
        vcf = sv_output_rule.output.vcf
    output:
        vcf = join(WORK_DIR, 'step10_{virus}_filter_vcf/breakpoints.vcf.gz'),
        tbi = join(WORK_DIR, 'step10_{virus}_filter_vcf/breakpoints.vcf.gz.tbi'),
    shell:
        "bcftools view {input.vcf} | "
        "egrep -e '^#|{wildcards.virus}' | "
        "bcftools filter -e \"IMPRECISE=1\" -Oz -o {output.vcf}"
        " && tabix -p vcf {output.vcf}"

# bcftools doesn't work with 4th columns in bed
# vcfanno doesn't work when mulitple genes overlaping one variant
# so using bedtools intersect plus a loop over the VCF
def overlap_with_genes(vcf_path, output_vcf_path, genes_path, work_dir, anno_name, anno_desc):
    overlap_tsv = join(work_dir, 'overlap.tsv')
    shell(f'bedtools intersect -a {vcf_path} -b {genes_path} -loj > {overlap_tsv}')
    gene_by_id = defaultdict(list)
    with open(overlap_tsv) as f:
        for l in f:
            fields = l.strip().split('\t')
            var_id = fields[2]
            VCF_COL_NUM = 10
            chrom, start, end, gene = fields[VCF_COL_NUM:VCF_COL_NUM+4]
            if gene.strip():
                gene_by_id[var_id].append(gene)

    ungz, gz = get_ungz_gz(output_vcf_path)
    with open_gzipsafe(vcf_path) as inp_f, open(ungz, 'w') as out_f:
        for l in inp_f:
            if l.startswith('##'):
                out_f.write(l)
            elif l.startswith('#'):
                out_f.write(f'##INFO=<ID={anno_name},Number=.,Type=String,Description="{anno_desc}">\n')
                out_f.write(l)
            else:
                fields = l.strip().split('\t')
                var_id = fields[2]
                info = fields[7]
                info_d = dict(kv.split('=') if '=' in kv else (kv, True) for kv in info.split(';'))
                mate_id = info_d.get('MATEID')
                genes = [g.replace('"', '').strip() for g in
                         gene_by_id.get(var_id, []) + gene_by_id.get(mate_id, [])
                         if g.replace('"', '').strip() != '.' and g is not None]
                if genes:
                    info += f';{anno_name}={",".join(genes)}'
                fields[7] = info
                out_f.write('\t'.join(fields) + '\n')
    shell(f'bgzip {ungz} && tabix -p vcf {gz}')

    before = count_vars(vcf_path)
    after = count_vars(output_vcf_path)
    assert before == after, (before, after)

rule annotate_with_viral_genes:
    input:
        vcf = rules.filter_vcf.output.vcf,
        tbi = rules.filter_vcf.output.tbi,
    output:
        vcf = join(WORK_DIR, 'step11_{virus}_viral_genes/breakpoints.viral_genes.vcf.gz'),
        tbi = join(WORK_DIR, 'step11_{virus}_viral_genes/breakpoints.viral_genes.vcf.gz.tbi'),
    params:
        work_dir = lambda wc, input, output: dirname(output.vcf)
    run:
        genes_bed = join(package_path(), 'data', f'{wildcards.virus}.bed')
        if not isfile(genes_bed):
            warn(f'No genes data for virus {wildcards.virus}, skipping annotation')
            shell(f'cp {input.vcf} {output.vcf}')
        else:
            overlap_with_genes(input.vcf, output.vcf, genes_bed, params.work_dir,
                'ViralGenes',
                'Viral genes that this breakpoint overlaps (and likely disrupts)')

host_genes_bed = join(package_path(), 'data', 'hg38_noalt.genes.bed.gz')
if GTF_PATH:  # user provided GTF file - overriding the default annotation BED
    rule prep_gtf:
        input:
            gtf_path = GTF_PATH,
            fai = rules.create_combined_reference.output.combined_fa + '.fai',
        output:
            gtf = join(WORK_DIR, 'host_genes_prep/hg38_noalt.genes.gtf'),
        params:
            work_dir = join(WORK_DIR, 'host_genes_prep'),
        run:
            # in order to pad genes by 100kb, we first need to subset the GTF to main chromosomes, and
            # remove chr prefixes from the fai file for bedtools.
            hg38_fai_to_grch38 = join(params.work_dir, 'grch38_noalt.fai')
            shell(f"cat {input.fai} | grep chr | sed 's/chrM/MT/' | sed 's/chr//' > {hg38_fai_to_grch38}")

            grch38_noalt_bed = join(params.work_dir, 'grch38_noalt.bed')
            shell(f"cat {hg38_fai_to_grch38} | awk '{{printf(\"%s\\t0\\t%d\\n\", $1, $2, $2-1)}}'"
                  f" > {grch38_noalt_bed}")

            shell(f"bedtools intersect -a {input.gtf_path} -b {grch38_noalt_bed}"
                  f" | grep -w gene | sed 's/^MT/M/' | sed 's/^/chr/' > {output.gtf}")

    rule gtf_to_bed:
        input:
            gtf = rules.prep_gtf.output.gtf,
        output:
            bed = join(WORK_DIR, 'host_genes_prep/hg38_noalt.genes.bed'),
        run:
            with open(input.gtf) as i_f, open(output.bed, 'w') as o_f:
                for l in i_f:
                    chrom, _, _, start, end, _, strand, _, attrs = l.strip().split('\t')
                    attr_d = dict(item.strip().split(' ') for item in attrs.split(';') if item)
                    gene = attr_d.get('gene_name').replace('"', '').strip()
                    if gene:
                        o_f.write(f'{chrom}\t{start}\t{end}\t{gene}\t.\t{strand}\n')

    host_genes_bed = rules.gtf_to_bed.output.bed

rule annotate_with_disrupted_host_genes:
    input:
        vcf = rules.annotate_with_viral_genes.output.vcf,
        bed = host_genes_bed,
        fai = rules.create_combined_reference.output.combined_fa + '.fai',
    output:
        vcf = join(WORK_DIR, 'step12_{virus}_host_genes_disrupted/breakpoints.annotated.vcf.gz'),
    params:
        work_dir = lambda wc, input, output: dirname(output.vcf)
    run:
        overlap_with_genes(input.vcf, output.vcf, input.bed, params.work_dir,
            'DisruptedGenes',
            'Host genes overlapping the breakpoint, that are probably disrupted '
            'by the viral integration event')

HOST_GENES_BASES_UPSTREAM = 100_000

rule slop_host_bed:
    input:
        bed = host_genes_bed,
        fai = rules.create_combined_reference.output.combined_fa + '.fai',
    output:
        bed = join(WORK_DIR, 'host_genes_prep/hg38_noalt_genes_slopped.bed'),
    params:
        work_dir = join(WORK_DIR, 'host_genes_prep'),
        bases_upstream = HOST_GENES_BASES_UPSTREAM,
    shell:
        # pad genes by 100kb
        # -l - The number of base pairs to subtract from the start coordinate
        # -s - Define -l and -r based on strand
        'zgrep -f <(cut -f1 {input.fai}) {input.bed} | '
        'bedtools slop -i stdin -l {params.bases_upstream} -r 0 -s -g {input.fai}'
        ' > {output.bed}'

rule annotate_with_host_genes_upstream:
    input:
        vcf = rules.annotate_with_disrupted_host_genes.output.vcf,
        bed = rules.slop_host_bed.output.bed,
        fai = rules.create_combined_reference.output.combined_fa + '.fai',
    output:
        vcf = join(WORK_DIR, 'step13_{virus}_host_genes_upstream/breakpoints.genes.host_cancer_genes.vcf.gz'),
    params:
        work_dir = lambda wc, input, output: dirname(output.vcf),
        bases_upstream = HOST_GENES_BASES_UPSTREAM,
    run:
        overlap_with_genes(input.vcf, output.vcf, input.bed, params.work_dir,
            f'GenesWithin{params.bases_upstream // 1000}kb',
            f'Genes that start within the {params.bases_upstream // 1000}kb distance '
            f'of the breakpoint. Based on https://www.biorxiv.org/content/10.1101/2020.02.12.942755v1, '
            f'that reported significant increases in chromatin accessibility exclusively '
            f'within 100 kbp of HPV integration sites.')

def merge_viruses_input_fn(wildcards):
    if not VIRUSES:
        selected_viruses_tsv = checkpoints.select_viruses.get(**wildcards).output.selected_viruses_tsv
        viruses = [v.strip() for v in open(selected_viruses_tsv).readlines() if v.strip()]
    else:
        viruses = VIRUSES
    return expand(rules.annotate_with_host_genes_upstream.output[0], virus=viruses)

rule merged_viruses:
    input:
        merge_viruses_input_fn
    output:
        RESULT_PATH
    run:
        if len(input) == 0:
            shell(f'touch {output}')
        else:
            if len(input) > 1:
                shell(f'bcftools merge {input} -Oz -o {output}')
            else:
                shell(f'cp {input} {output}')
            shell(f'tabix -p vcf {output}')
