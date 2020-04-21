from collections import defaultdict
from os.path import isfile, join, basename, dirname, abspath
import subprocess
from ngs_utils.file_utils import which, safe_mkdir, verify_file, get_ungz_gz
from ngs_utils.logger import warn
from hpc_utils import hpc
from ngs_utils.vcf_utils import add_cyvcf2_hdr
from ngs_utils.reference_data import get_key_genes_bed
from ngs_utils.vcf_utils import count_vars, vcf_contains_field, iter_vcf

from oncoviruses import package_path


#############################
#### Reading parameters #####

INPUT_BAM = config['input_bam']
OUTPUT_DIR = abspath(config['output_dir'])
RESULT_PATH = join(OUTPUT_DIR, 'breakpoints.vcf.gz')
PRIO_TSV = join(OUTPUT_DIR, 'prioritized_oncoviruses.tsv')


SAMPLE = config['sample_name']
VIRUS = config.get('virus', None)
if VIRUS:
    VIRUS = VIRUS.upper()
THREADS = config.get('cores', 10)

if config.get('genomes_dir'):
    hpc.set_genomes_dir(config.get('genomes_dir'))

GENOME = 'hg38'
COMBINED_FA = config.get('combined_fa') or \
    hpc.get_ref_file(genome=GENOME, key='fa_plus_gdc_viruses', must_exist=False)
HOST_FA = config.get('host_fa') or \
    hpc.get_ref_file(genome=GENOME, key='fa', must_exist=False)
VIRUSES_FA = config.get('viruses_fa') or \
    hpc.get_ref_file(key='gdc_viral_fa', must_exist=False)
assert VIRUSES_FA and (COMBINED_FA or HOST_FA)

WORK_DIR = join(OUTPUT_DIR, 'work')

SV_CALLER = 'manta'

PY2_ENV_PATH = config.get('py2_env_path')
py2_conda_cmd = ''
if PY2_ENV_PATH:
    py2_conda_cmd = 'export PATH=' + PY2_ENV_PATH + '/bin:$PATH; '


rule all:
    input:
        RESULT_PATH
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
        params:
            rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
        shell:
            "bwa mem -Y -t{threads} -R '{params.rg}' {input.gdc_fa} {input.fq1} {input.fq2} "
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

    ONCOVIRAL_SOURCE_URL = 'https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files'
    rule prioritize_viruses:
        input:
            mosdepth_regions_bed_gz = rules.mosdepth.output.mosdepth_regions_bed_gz,
            mosdepth_thresholds_bed_gz = rules.mosdepth.output.mosdepth_thresholds_bed_gz,
        output:
            oncoviruses_tsv = PRIO_TSV,
        shell:
            "echo '## Viral sequences (from {ONCOVIRAL_SOURCE_URL}) found in unmapped reads' > {output.oncoviruses_tsv} &&"
            "echo '## Sample: {SAMPLE}' >> {output.oncoviruses_tsv} && "
            "echo '#virus\tsize\tdepth\t1x\t5x\t25x' >> {output.oncoviruses_tsv} && "
            "paste <(gunzip -c {input.mosdepth_regions_bed_gz}) <(zgrep -v ^# {input.mosdepth_thresholds_bed_gz}) | "
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
                    f'Not found any viruses with significant coverage. '
                    f'You can explore the full list at {input.tsv} and rerun with '
                    f'the --virus option explicitly.')
                shell('touch {output.selected_viruses_tsv}')
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
    params:
        rg = f'@RG\\tID:{SAMPLE}\\tSM:{SAMPLE}'
    shell:
        # using the polyidus bwa command.
        # -T1 = minimum score to output [default 30]
        # -a  = output all alignments for SE or unpaired PE
        # -Y  = use soft clipping for supplementary alignments
        "bwa mem -a -Y -t{threads} -R '{params.rg}' {input.virus_bwa_prefix} {input.fq1} {input.fq2}"
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

# if SV_CALLER == 'manta':
MANTA_IMG = 'quay.io/biocontainers/manta:1.6.0--py27_0'
rule run_manta:
    input:
        bam = rules.bwa_viral_bridging_to_comb_ref.output.comb_bam_possorted,
        bai = rules.index_comb_ref_bam.output.bai,
        ref = COMBINED_FA if COMBINED_FA is not None else rules.create_combined_reference.output.combined_fa,
    output:
        vcf = join(WORK_DIR, 'step8_{virus}_manta/results/variants/candidateSV.vcf.gz'),
    params:
        work_dir = join(WORK_DIR, 'step8_{virus}_manta'),
        image = MANTA_IMG,
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
            f'{run_script} -m local'
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
            shell(f'{py2_conda_cmd} {tool_cmd}')

# rule run_manta:
#     input:
#         run_script = rules.prep_manta.output.run_script,
#     output:
#         vcf = join(WORK_DIR, 'step8_{virus}_manta/results/variants/candidateSV.vcf.gz'),
#     params:
#         image = MANTA_IMG,
#     threads: max(10, THREADS)
#     group: 'manta'
#     run:
#         tool_cmd = f'{input.run_script} -m local'
#
#         if subprocess.run(f'docker images -q {params.image} 2>/dev/null').returncode == 0:
#             bam_dir = abspath(dirname(input.bam))
#             ref_dir = abspath(dirname(input.ref))
#             shell(
#                 f'docker run '
#                 f'-v{bam_dir}:{bam_dir} '
#                 f'-v{ref_dir}:{ref_dir} '
#                 f'-v{params.work_dir}:/work '
#                 f'{params.image} {tool_cmd.replace(params.work_dir, "/work")}'
#             )
#         else:
#             shell(f'{py2_conda_cmd} {tool_cmd}')


sv_output_rule = rules.run_manta

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

# rule snpeff:
#     input:
#         vcf = rules.filter_vcf.output.vcf,
#     output:
#         vcf = join(WORK_DIR, 'step10_{virus}_snpeff/breakpoints.eff.vcf.gz')
#     params:
#         genome = GENOME,
#         tmp_dir = join(WORK_DIR, 'step10_{virus}_snpeff', 'tmp'),
#         csv  = join(WORK_DIR, 'step10_{virus}_snpeff', 'snpeff-stats.csv'),
#         html = join(WORK_DIR, 'step10_{virus}_snpeff', 'snpeff-stats.html'),
#     resources:
#         mem_mb = 2000
#     run:
#         safe_mkdir(params.tmp_dir)
#         snpeff_db = hpc.get_ref_file(genome=params.genome, key='snpeff')
#         snpeff_db_dir = dirname(snpeff_db)
#         assert params.genome == 'hg38', 'Only hg38 is suported for snpEff'
#         snpeff_db_name = 'GRCh38.86'  # for some reason it doesn't matter if the subdir is named GRCh38.92
#
#         jvm_opts = f'-Xms750m -Xmx4g'
#         java_args = f'-Djava.io.tmpdir={params.tmp_dir}'
#
#         shell('snpEff {jvm_opts} {java_args} '
#               '-dataDir {snpeff_db_dir} {snpeff_db_name} '
#               '-hgvs -cancer -i vcf -o vcf '
#               '-csvStats {params.csv} -s {params.html} '
#               '{input.vcf} '
#               '| bgzip -c > {output.vcf} && tabix -p vcf {output.vcf}')
#         verify_file(output.vcf, is_critical=True)

# bcftools doesn't work with 4th columns in bed
# vcfanno doesn't work when mulitple genes overlaping one variant
# so using bedtools intersect plus a cyvcf2 loop
rule annotate_with_viral_genes:
    input:
        vcf = rules.filter_vcf.output.vcf,
        tbi = rules.filter_vcf.output.tbi,
    output:
        vcf = join(WORK_DIR, 'step10_{virus}_viral_genes/breakpoints.viral_genes.vcf.gz'),
        tbi = join(WORK_DIR, 'step10_{virus}_viral_genes/breakpoints.viral_genes.vcf.gz.tbi'),
    params:
        overlap_tsv = join(WORK_DIR, 'step10_{virus}_viral_genes/overlap.tsv')
    run:
        genes_bed = join(package_path(), 'data', f'{wildcards.virus}.bed')
        if not isfile(genes_bed):
            warn(f'No genes data for virus {wildcards.virus}, skipping annotation')
            shell('cp {input.vcf} {output.vcf}')
        else:
            shell('bedtools intersect -a {input.vcf} -b {genes_bed} -loj > {params.overlap_tsv}')

            gene_by_id = defaultdict(list)
            with open(params.overlap_tsv) as f:
                for l in f:
                    _, _, var_id, _, _, _, _, _, _, _, _, gene = l.strip().split('\t')
                    gene_by_id[var_id].append(gene)

            def proc_rec(rec, vcf):
                genes = [g for g in
                         gene_by_id.get(rec.ID, []) + gene_by_id.get(rec.INFO.get('MATEID', '.'), [])
                         if g != '.']
                if genes:
                    rec.INFO['ViralGenes'] = ','.join(genes)
                return rec
            def proc_hdr(vcf):
                add_cyvcf2_hdr(vcf, 'ViralGenes', '.', 'String',
                               'Viral genes that this breakpoint overlaps (and likely disrupts)')
            iter_vcf(input.vcf, output.vcf, proc_rec=proc_rec, proc_hdr=proc_hdr)

        before = count_vars(input.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)

rule annotate_with_host_genes:
    input:
        vcf = rules.annotate_with_viral_genes.output.vcf,
        key_genes_bed = get_key_genes_bed(GENOME),
        # coding_regions_bed = hpc.get_ref_file(GENOME, key='coding_regions'),
        fai = COMBINED_FA + '.fai',
    output:
        vcf = join(WORK_DIR, 'step11_{virus}_host_genes/breakpoints.genes.host_cancer_genes.vcf.gz'),
    params:
        overlap_tsv = join(WORK_DIR, 'step11_{virus}_host_genes/overlap.tsv')
    run:
        # -D b Report the closest feature in B, and its distance to A as an extra column. Use negative distances
        #      to report upstream features. -Db means report distance with respect to B: When B is on the - strand,
        #      “upstream” means A has a higher (start,stop).
        shell('bedtools closest -D b -a {input.vcf} -b {input.key_genes_bed} -g {input.fai}'
              ' > {params.overlap_tsv}')
        gene_by_id = defaultdict(list)
        with open(params.overlap_tsv) as f:
            for l in f:
                try:
                    _, _, var_id, _, _, _, _, _, _, _, _, gene, _, strand, distance = l.strip().split('\t')
                    distance = int(distance)
                except:
                    pass
                else:
                    # Reporting either overlaps (distance ~ 0), or upstream events that might
                    # affect gene expression:
                    if distance < 50 or strand == '+' and distance < 0 or strand == '-' and distance > 0:
                        gene_by_id[var_id].append(gene)

        def proc_rec(rec, vcf):
            genes = [g for g in
                     gene_by_id.get(rec.ID, []) + gene_by_id.get(rec.INFO.get('MATEID', '.'), [])
                     if g != '.']
            if genes:
                rec.INFO['NearestCancerGenes'] = ','.join(genes)
            return rec
        def proc_hdr(vcf):
            add_cyvcf2_hdr(vcf, 'NearestCancerGenes', '.', 'String',
                           'Nearest cancer gene to the breakpoint')
        iter_vcf(input.vcf, output.vcf, proc_rec=proc_rec, proc_hdr=proc_hdr)

        before = count_vars(input.vcf)
        after = count_vars(output.vcf)
        assert before == after, (before, after)

last_vcf_rule = rules.annotate_with_host_genes

# rule tabix_result_vcf:
#     input:
#         vcf = rules.annotate_with_host_genes.output.vcf
#     output:
#         tbi = rules.annotate_with_host_genes.output.vcf + '.tbi'
#     shell:
#         'tabix -p vcf {input}'

def merge_viruses_input_fn(wildcards):
    if not VIRUS:
        selected_viruses_tsv = checkpoints.select_viruses.get(**wildcards).output.selected_viruses_tsv
        viruses = [v.strip() for v in open(selected_viruses_tsv).readlines() if v.strip()]
    else:
        viruses = [VIRUS]
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


# TODO: count viruses, breakpoints and genes for MultiQC report
