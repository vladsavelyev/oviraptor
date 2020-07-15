import os
from collections import defaultdict
from os.path import isfile, join, basename, dirname, abspath
import subprocess
import glob
import re
from ngs_utils.file_utils import which, safe_mkdir, verify_file, get_ungz_gz
from ngs_utils.logger import warn, critical
from ngs_utils.vcf_utils import add_cyvcf2_hdr
from ngs_utils.reference_data import get_key_genes_bed
from ngs_utils.vcf_utils import count_vars, vcf_contains_field, iter_vcf
from reference_data import api as refdata
from oncoviruses import package_path


#############################
#### Reading parameters #####

INPUT_BAM = config['input_bam']
OUTPUT_DIR = abspath(config['output_dir'])
RESULT_PATH = join(OUTPUT_DIR, 'breakpoints.vcf.gz')
PRIO_TSV = join(OUTPUT_DIR, 'prioritized_oncoviruses.tsv')


SAMPLE = config['sample_name']
VIRUSES = config.get('viruses', None)
if VIRUSES:
    VIRUSES = [v.upper() for v in VIRUSES.split(',')]
ONLY_DETECT = config.get('only_detect', False)

THREADS = config.get('cores', 10)

if config.get('genomes_dir'):
    refdata.find_genomes_dir(config.get('genomes_dir'))

GENOME = 'hg38'
COMBINED_FA = config.get('combined_fa') or \
    refdata.get_ref_file(genome=GENOME, key='fa_plus_gdc_viruses', must_exist=False)
HOST_FA = config.get('host_fa') or \
    refdata.get_ref_file(genome=GENOME, key='fa', must_exist=False)
VIRUSES_FA = config.get('viruses_fa') or \
    refdata.get_ref_file(genome=GENOME, key='gdc_viral_fa', must_exist=False)
assert VIRUSES_FA and (COMBINED_FA or HOST_FA)

WORK_DIR = join(OUTPUT_DIR, 'work')

SV_CALLER = 'lumpy'

GTF_PATH = config.get('gtf_file')
if not GTF_PATH:
    pyens = refdata.get_ref_file(GENOME, key='pyensembl_data', must_exist=False)
    if pyens:
        pattern = join(refdata.get_ref_file(GENOME, key='pyensembl_data', must_exist=False),
                               'GRCh38/ensembl*/Homo_sapiens.GRCh38.*.gtf.gz')
        found = glob.glob(pattern)
        if not found:
            critical(f'Could not find a GTF file in pyensembl: {pattern}')
        GTF_PATH = found[-1]
    else:
        warn('GTF file not provided, skipping annotating with host genes. Use --gtf option if you want annotations.')


rule all:
    input:
        RESULT_PATH if not ONLY_DETECT else [],
        PRIO_TSV if ONLY_DETECT or not VIRUSES else [],
    # input: join(WORK_DIR, 'all.done')


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
            shell("bwa mem -Y -t{threads} -R '{params.rg}' {input.gdc_fa} {input.fq1} {input.fq2} "
                  " | samtools sort -@{threads} -Obam -o {output.gdc_bam}")

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
            image = 'quay.io/biocontainers/mosdepth:0.2.9--hbeb723e_0'
        threads: THREADS
        run:
            chroms_bed = join(params.work_dir, "chroms.bed")
            awk_cmd = "awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}'"
            shell(f'{awk_cmd} {input.fai} > {chroms_bed}')
            mosdepth_cmd = (
                f'mosdepth {params.prefix} {input.bam} -t{threads} -n --thresholds 1,5,25 --by {chroms_bed} '
            )
            if subprocess.run(f'docker images -q {params.image} 2>/dev/null', shell=True).returncode == 0:
                bam_dir = abspath(dirname(input.bam))
                shell(
                    f'docker run -v{bam_dir}:{bam_dir} -v{params.work_dir}:/work '
                    f'{params.image} ' +
                    f'{mosdepth_cmd.replace(params.work_dir, "/work")}'
                )
            else:
                shell(mosdepth_cmd)

    # we need at least one of these conditions to call significance:
    MIN_1x_PCT = 50.0  #  % of the viral sequence that must be covered at at least 1x (probably non-integrating)
    MIN_5x_LEN = 300   #  viral base pairs must be covered at at least 5x (which is amplified, thus integrating)
    ONCOVIRAL_SOURCE_URL = 'https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files'
    rule prioritize_viruses:
        input:
            mosdepth_regions_bed_gz = rules.mosdepth.output.mosdepth_regions_bed_gz,
            mosdepth_thresholds_bed_gz = rules.mosdepth.output.mosdepth_thresholds_bed_gz,
        output:
            oncoviruses_tsv = PRIO_TSV,
        params:
            completeness_share = MIN_1x_PCT / 100.0
        shell:
            "echo '## Viral sequences (from {ONCOVIRAL_SOURCE_URL}) found in unmapped reads' > {output.oncoviruses_tsv} &&"
            "echo '## Sample: {SAMPLE}' >> {output.oncoviruses_tsv} && "
            "echo '## Minimal completeness: {MIN_1x_PCT}% at 1x or {MIN_5x_LEN}bp at 5x' >> {output.oncoviruses_tsv} && "
            "echo '#virus\tsize\tdepth\t1x\t5x\t25x\tsignificance' >> {output.oncoviruses_tsv} && "
            "paste <(gunzip -c {input.mosdepth_regions_bed_gz}) <(zgrep -v ^# {input.mosdepth_thresholds_bed_gz}) | "
            "awk 'BEGIN {{FS=\"\\t\"}} {{ printf(\"%s\\t%d\\t%3.1f\\t%3.3f\\t%3.3f\\t%3.3f\\t%s\\n\", "
            "$1, $3, $4, $9/$3, $10/$3, $11/$3, (($10>{MIN_5x_LEN} || $11/$3>{params.completeness_share}) ? \"significant\" : \".\")) }}' | "
            "sort -n -r -k5,5 -k6,6 -k4,4 -k3,3 >> {output.oncoviruses_tsv}"

    checkpoint select_viruses:
        input:
            tsv = rules.prioritize_viruses.output.oncoviruses_tsv,
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
                shell('touch {output.selected_viruses_tsv}')
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
        virus_bwt = join(WORK_DIR, '{virus}', '{virus}.fa.bwt'),  # one of bwa index files
    shell:
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
    params:
        rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
    shell:
        # using the polyidus bwa command.
        # -T1 = minimum score to output [default 30]
        # -a  = output all alignments for SE or unpaired PE
        # -Y  = use soft clipping for supplementary alignments
        "bwa mem -a -Y -t{threads} -R '{params.rg}' {input.virus_bwa_prefix} {input.fq1} {input.fq2}"
        " | samtools sort -@{threads} -Obam -o {output.virus_bam_possorted}"


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
    params:
        rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
    shell:
        # using the polyidus bwa command.
        # -T1 = minimum score to output [default 30]
        # -a  = output all alignments for SE or unpaired PE
        # -Y  = use soft clipping for supplementary alignments
        "bwa mem -a -Y -t{threads} -R '{params.rg}' {input.bwa_prefix} {input.fq1} {input.fq2}"
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
        fai = (COMBINED_FA if COMBINED_FA is not None else rules.create_combined_reference.output.combined_fa) + '.fai',
    output:
        vcf = join(WORK_DIR, 'step8_{virus}_lumpy.vcf'),
    params:
        lumpy = join(package_path(), 'lumpy', 'lumpy'),
    group: 'lumpy'
    run:
        if refdata.ref_file_exists(GENOME, 'blacklist'):
            blacklist = refdata.get_ref_file(GENOME, 'blacklist')
            blacklist_opt = f'{f"-x {blacklist} " if isfile(blacklist) else ""}'
        shell(
            '{params.lumpy} '
            '-mw 4 '
            '-tt 0 '
            '-pe id:sample,bam_file:{input.disc},histo_file:{input.histo},mean:500,stdev:50,read_length:151,min_non_overlap:151,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 '
            '-sr id:sample,bam_file:{input.split},back_distance:10,weight:1,min_mapping_threshold:20 '
            '| gsort - {input.fai} > {output.vcf}'
        )

rule run_manta:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
        bai = rules.index_comb_ref_bam.output.bai,
        ref = COMBINED_FA if COMBINED_FA is not None else rules.create_combined_reference.output.combined_fa,
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
    sv_output_rule = rules.run_lumpy

rule filter_vcf:
    input:
        vcf = sv_output_rule.output.vcf
    output:
        vcf = join(WORK_DIR, 'step9_{virus}_manta_filter/breakpoints.vcf.gz'),
        tbi = join(WORK_DIR, 'step9_{virus}_manta_filter/breakpoints.vcf.gz.tbi'),
    shell:
        "bcftools view {input.vcf} | "
        "egrep -e '^#|{wildcards.virus}' | "
        "bcftools filter -e \"IMPRECISE=1\" -Oz -o {output.vcf}"
        " && tabix -p vcf {output.vcf}"

# bcftools doesn't work with 4th columns in bed
# vcfanno doesn't work when mulitple genes overlaping one variant
# so using bedtools intersect plus a cyvcf2 loop
def overlap_with_genes(vcf_path, output_vcf_path, genes_path, work_dir, anno_name, anno_desc):
    overlap_tsv = join(work_dir, 'overlap.tsv')
    shell(f'bedtools intersect -a {vcf_path} -b {genes_path} -loj > {overlap_tsv}')
    gene_by_id = defaultdict(list)
    with open(overlap_tsv) as f:
        for l in f:
            fields = l.strip().split('\t')
            var_id = fields[2]
            VCF_COL_NUM = 10
            gene = None
            if genes_path.endswith('.bed') or genes_path.endswith('.bed.gz'):
                chrom, start, end, gene = fields[VCF_COL_NUM:VCF_COL_NUM+4]
            else:  # GTF
                chrom, _, _, start, end, _, strand, _, info = fields[VCF_COL_NUM:]
                m = re.match(r'.*gene_name "(\S+)";.*', info)
                if m:
                    gene = m.group(1)
            if gene:
                gene_by_id[var_id].append(gene)

    def proc_rec(rec, vcf):
        genes = [g for g in
                 gene_by_id.get(rec.ID, []) + gene_by_id.get(rec.INFO.get('MATEID', '.'), [])
                 if g != '.' and g is not None]
        if genes:
            rec.INFO[anno_name] = ','.join(genes)
        return rec
    def proc_hdr(vcf):
        add_cyvcf2_hdr(vcf, anno_name, '.', 'String', anno_desc)
    iter_vcf(vcf_path, output_vcf_path, proc_rec=proc_rec, proc_hdr=proc_hdr)

    before = count_vars(vcf_path)
    after = count_vars(output_vcf_path)
    assert before == after, (before, after)

rule annotate_with_viral_genes:
    input:
        vcf = rules.filter_vcf.output.vcf,
        tbi = rules.filter_vcf.output.tbi,
    output:
        vcf = join(WORK_DIR, 'step10_{virus}_viral_genes/breakpoints.viral_genes.vcf.gz'),
        tbi = join(WORK_DIR, 'step10_{virus}_viral_genes/breakpoints.viral_genes.vcf.gz.tbi'),
    params:
        work_dir = join(WORK_DIR, 'step10_{virus}_viral_genes')
    run:
        genes_bed = join(package_path(), 'data', f'{wildcards.virus}.bed')
        if not isfile(genes_bed):
            warn(f'No genes data for virus {wildcards.virus}, skipping annotation')
            shell('cp {input.vcf} {output.vcf}')
        else:
            overlap_with_genes(input.vcf, output.vcf, genes_bed, params.work_dir,
                'ViralGenes',
                'Viral genes that this breakpoint overlaps (and likely disrupts)')

last_vcf_rule = rules.annotate_with_viral_genes

if GTF_PATH:
    rule prep_gtf:
        input:
            gtf_path = GTF_PATH,
            fai = COMBINED_FA + '.fai',
        output:
            gtf = join(WORK_DIR, 'step11_{virus}_host_genes/hg38_noalt.genes.gtf'),
        params:
            work_dir = join(WORK_DIR, 'step11_{virus}_host_genes'),
        run:
            # in order to pad genes by 100kb, we first need to subset the GTF to main chromosomes, and
            # remove chr prefixes from the fai file for bedtools.
            hg38_fai_to_grch38 = join(params.work_dir, 'grch38_noalt.fai')
            shell("cat {input.fai} | grep chr | sed 's/chrM/MT/' | sed 's/chr//' > {hg38_fai_to_grch38}")

            grch38_noalt_bed = join(params.work_dir, 'grch38_noalt.bed')
            shell("cat {hg38_fai_to_grch38} | awk '{{printf(\"%s\\t0\\t%d\\n\", $1, $2, $2-1)}}' > {grch38_noalt_bed}")

            shell("bedtools intersect -a {input.gtf_path} -b {grch38_noalt_bed}"
                  " | grep -w gene | sed 's/^MT/M/' | sed 's/^/chr/' > {output.gtf}")

    rule annotate_with_disrupted_host_genes:
        input:
            vcf = rules.annotate_with_viral_genes.output.vcf,
            gtf = rules.prep_gtf.output.gtf,
            fai = COMBINED_FA + '.fai',
        output:
            vcf = join(WORK_DIR, 'step11_{virus}_host_genes/breakpoints.annotated.vcf.gz'),
        params:
            work_dir = join(WORK_DIR, 'step11_{virus}_host_genes'),
        run:
            overlap_with_genes(input.vcf, output.vcf, input.gtf, params.work_dir,
                'DisruptedGenes',
                'Host genes overlapping the breakpoint, that are probably disrupted '
                'by the viral integration event')

    HOST_GENES_BASES_UPSTREAM = 100_000

    rule slop_gtf:
        input:
            gtf = rules.prep_gtf.output.gtf,
            fai = COMBINED_FA + '.fai',
        output:
            gtf = join(WORK_DIR, 'step12_{virus}_host_genes/hg38_noalt.genes.slopped.gtf'),
        params:
            work_dir = join(WORK_DIR, 'step12_{virus}_host_genes'),
            bases_upstream = HOST_GENES_BASES_UPSTREAM,
        shell:
            # pad genes by 100kb
            # -l - The number of base pairs to subtract from the start coordinate
            # -s - Define -l and -r based on strand
            'bedtools slop -i {input.gtf} -l {params.bases_upstream} -r 0 -s -g {input.fai}'
            ' > {output.gtf}'

    rule annotate_with_host_genes_upstream:
        input:
            vcf = rules.annotate_with_viral_genes.output.vcf,
            gtf = rules.slop_gtf.output.gtf,
            fai = COMBINED_FA + '.fai',
        output:
            vcf = join(WORK_DIR, 'step12_{virus}_host_genes/breakpoints.genes.host_cancer_genes.vcf.gz'),
        params:
            work_dir = join(WORK_DIR, 'step12_{virus}_host_genes'),
            bases_upstream = HOST_GENES_BASES_UPSTREAM,
        run:
            overlap_with_genes(input.vcf, output.vcf, input.gtf, params.work_dir,
                f'GenesWithin{params.bases_upstream // 1000}kb',
                f'Genes that start within the {params.bases_upstream // 1000}kb distance '
                f'of the breakpoint. Based on https://www.biorxiv.org/content/10.1101/2020.02.12.942755v1, '
                f'that reported significant increases in chromatin accessibility exclusively '
                f'within 100 kbp of HPV integration sites.')

    last_vcf_rule = rules.annotate_with_host_genes_upstream

def merge_viruses_input_fn(wildcards):
    if not VIRUSES:
        selected_viruses_tsv = checkpoints.select_viruses.get(**wildcards).output.selected_viruses_tsv
        viruses = [v.strip() for v in open(selected_viruses_tsv).readlines() if v.strip()]
    else:
        viruses = VIRUSES
    return expand(last_vcf_rule.output[0], virus=viruses)

rule merged_viruses:
    input:
        merge_viruses_input_fn
    output:
        RESULT_PATH
    run:
        if len(input) == 0:
            shell('touch {output}')
        else:
            if len(input) > 1:
                shell('bcftools merge {input} -Oz -o {output}')
            else:
                shell('cp {input} {output}')
            shell('tabix -p vcf {output}')
