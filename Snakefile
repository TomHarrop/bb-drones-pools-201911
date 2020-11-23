#!/usr/bin/env python3

import multiprocessing
import pandas


def get_min_cutoff(wildcards):
    cutoff_file = checkpoints.genotype.get(**wildcards).output['cutoffs']
    cutoffs = pandas.read_csv(cutoff_file,
                              header=None,
                              index_col=0)
    return cutoffs.loc['min_depth', 1]


def get_max_cutoff(wildcards):
    cutoff_file = checkpoints.genotype.get(**wildcards).output['cutoffs']
    cutoffs = pandas.read_csv(cutoff_file,
                              header=None,
                              index_col=0)
    return cutoffs.loc['max_depth', 1]

###########
# GLOBALS #
###########

bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.6')
r = ('shub://TomHarrop/r-containers:r_3.6.1'
     '@e1eb426cd153fd0669bc24508673228d2f25dd76')
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
vcftools = ('shub://TomHarrop/variant-utils:vcftools_0.1.16'
            '@d64cc5a37951760be575c43024c66e69b2563166')
whatshap = 'shub://TomHarrop/variant-utils:whatshap_0.18'
minimap = 'shub://TomHarrop/align-utils:minimap2_2.17r941'

ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
fai = f'{ref}.fai'

sample_info = 'data/combined_sampleinfo.csv'
cnv_map = 'data/cnv_map.txt'

########
# MAIN #
########

# read sample data
sample_df = pandas.read_csv(sample_info,
                            index_col='sample')

all_samples = sorted(set(sample_df.index))
pool_samples = [x for x in all_samples if x.endswith('_pool')]
drone_samples = [x for x in all_samples if x.endswith('_drone')]

drone_indivs = [x.split('_')[0] for x in drone_samples]
pool_indivs = [x.split('_')[0] for x in pool_samples]

both_indivs = [x for x in drone_indivs if x in pool_indivs]

# read reference data
fai_pd = pandas.read_csv(fai, sep='\t', header=None)
all_chr = sorted(set(fai_pd[0]))
autosomes = [x for x in all_chr if x.startswith('NC_')]

# map barcodes to names
# pool_data = pandas.read_csv(pools_csv,
#                             index_col='sample')
# pool_samples = sorted(set(pool_data.index))
# drone_data = pandas.read_csv(drones_csv,
#                              index_col='sample')
# drone_samples = sorted(set(drone_data.index))
# all_samples = {'drones': drone_samples,
#                'pools': pool_samples}

# # generate output files for demultiplexing
# sample_data = {'drones': drone_data,
#                'pools': pool_data}

# sample_df = pandas.concat(sample_data)

#########
# RULES #
#########

wildcard_constraints:
    chr = '|'.join(autosomes)

# rules
rule target:
    input:
        'output/020_filtered-genotypes/filtered.vcf.gz',

# rule phase_reads:
#     input:
#         vcf = 'output/030_phased/phased_pools.vcf.gz',
#         ref = 'output/010_genotypes/pools/015_ref/ref.fasta',
#         bam = 'output/010_genotypes/pools/merged.bam'
#     output:
#         bam = 'output/030_phased/phased_pools.bam'
#     log:
#         'output/logs/phase_reads.log'
#     singularity:
#         whatshap
#     shell:
#         'whatshap haplotag '
#         '-o {output.bam} '
#         '--reference {input.ref} '
#         '{input.vcf} '
#         '{input.bam} '
#         '&> {log} '
#         '|| true '  # don't know why it's failing but i need to see the output

# rule gtf_blocks:
#     input:
#         vcf = 'output/030_phased/phased_pools.vcf.gz',
#         fai = 'output/010_genotypes/pools/015_ref/ref.fasta.fai'
#     output:
#         gtf = 'output/040_phased-indivs/{indiv}.gtf',
#         stats = 'output/040_phased-indivs/{indiv}.txt',
#     log:
#         'output/logs/040_phased-indivs/{indiv}_gtf-blocks.log'
#     singularity:
#         whatshap
#     shell:
#         'whatshap stats '
#         '--gtf {output.gtf} '
#         '--chr-lengths <(cut -f1,2 {input.fai}) '
#         '--sample {wildcards.indiv} '
#         '{input.vcf} '
#         '> {output.stats} '
#         '2> {log}'


# rule phase_pools:
#     input:
#         vcf = 'output/020_filtered-genotypes/pools/filtered.vcf.gz',
#         bam = 'output/010_genotypes/pools/merged.bam',
#         ref = 'output/010_genotypes/pools/015_ref/ref.fasta',
#         fai = 'output/010_genotypes/pools/015_ref/ref.fasta.fai',
#         drone_vcf = 'output/030_phased/phased_drones.vcf.gz'
#     output:
#         'output/030_phased/phased_pools.vcf'
#     log:
#         'output/logs/phase_pools.log'
#     singularity:
#         whatshap
#     shell:
#         'whatshap phase '
#         '--reference {input.ref} '
#         '-o {output} '
#         '--indels '
#         '{input.vcf} '
#         '{input.bam} '
#         '{input.drone_vcf} ' # add phased drone vcf here
#         '&> {log}'

# # manually add a dummy phase set block to all drone calls
# rule add_phase_set:
#     input:
#         vcf = 'output/020_filtered-genotypes/drones/filtered.vcf'
#     output:
#         vcf = 'output/030_phased/phased_drones.vcf',
#         tmp = temp('output/030_phased/phased_drones_noheader.vcf')
#     log:
#         'output/logs/add_phase_set.log'
#     singularity:
#         r
#     shell:
#         'Rscript src/manual_vcf_parse.R {input.vcf} {output.tmp} '
#         '&> {log} ; '
#         'cat <(grep "^##" {input.vcf}) '
#         '<(echo \'##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\') '
#         '<(grep "^#[^#]" {input.vcf}) '
#         ' {output.tmp} '
#         ' > {output.vcf} '


rule split_target:
    input:
        expand('output/020_filtered-genotypes/{set}.renamed.vcf.gz',
               set=['pool', 'drone'])

rule indiv_target:
    input:
        expand('output/000_tmp/drones/{indiv}/{chr}.bam',
               indiv=['BB34'],
               chr=['NC_037640.1']),
        expand('output/000_tmp/pools/{indiv}/{chr}.bam',
               indiv=['BB34'],
               chr=['NC_037640.1'])

rule phase_target:
    input:
        expand('output/040_phased-indivs/{indiv}/{chr}.vcf',
               indiv=both_indivs,
               chr=autosomes)

# phase?
rule phase:
    input:
        pool_vcf = 'output/000_tmp/pools/{indiv}/{chr}.vcf.gz',
        pool_bam = 'output/000_tmp/pools/{indiv}/{chr}.bam',
        pool_bai = 'output/000_tmp/pools/{indiv}/{chr}.bam.bai',
        drone_bam = 'output/000_tmp/drones/{indiv}/{chr}.bam',
        drone_bai = 'output/000_tmp/drones/{indiv}/{chr}.bam.bai'
    output:
        'output/040_phased-indivs/{indiv}/{chr}.vcf'
    log:
        'output/logs/phase.{indiv}.{chr}.log'
    container:
        whatshap
    shell:
        'whatshap phase '
        '-o {output} '
        '{input.pool_vcf} '
        '{input.drone_bam} '
        '{input.pool_bam} '
        '&>{log}'


# rename the pool reads
rule rename_bam:
    input:
        'output/000_tmp/pools/{indiv}/{chr}.sam'
    output:
        'output/000_tmp/pools/{indiv}/{chr}.bam'
    container:
        samtools
    shell:
        'samtools addreplacerg '
        '-r "ID:{wildcards.indiv}" '
        '-r "SM:{wildcards.indiv}" '
        '-O BAM '
        '{input} '
        '> {output}'

rule extract_pool_bam:
    input:
        bam = 'output/010_genotypes/merged.bam'
    output:
        pipe('output/000_tmp/pools/{indiv}/{chr}.sam')
    log:
        'output/logs/extract_pool_bam.{indiv}.{chr}.log'
    params:
        query = lambda wildcards: f'{wildcards.indiv}_pool',
    container:
        samtools
    shell:
        'samtools view -h '
        '-r {params.query} '
        '{input.bam} '
        '{wildcards.chr} '
        '>>{output} '
        '2>{log}'

# map
rule sam_to_bam:
    input:
        'output/000_tmp/drones/{indiv}/{chr}.sam'
    output:
        bam = 'output/000_tmp/drones/{indiv}/{chr}.bam',
    container:
        samtools
    shell:
        'samtools view -bh {input} '
        '| samtools sort '
        '> {output.bam} '

rule map_read:
    input:
        ref = ref,
        read = 'output/000_tmp/drones/{indiv}/{chr}.renamed.fa'
    output:
        pipe('output/000_tmp/drones/{indiv}/{chr}.sam')
    log:
        'output/logs/map_read.{indiv}.{chr}.log'
    params:
        rg = lambda wildcards: f'@RG\\tID:{wildcards.indiv}\\tSM:{wildcards.indiv}'
    container:
        minimap
    shell:
        'minimap2 '
        '-ax asm5 '
        '-R \'{params.rg}\' '
        '{input.ref} '
        '{input.read} '
        '>>{output} '
        '2>{log}'


# change the name of the read to the name of the drone
rule rename_read:
    input:
        'output/000_tmp/drones/{indiv}/{chr}.fa'
    output:
        'output/000_tmp/drones/{indiv}/{chr}.renamed.fa'
    params:
        query = lambda wildcards: f's/>{wildcards.chr}/>{wildcards.indiv}/g'
    container:
        samtools
    shell:
        'sed '
        '\'{params.query}\' '
        '{input} '
        '>{output}'

# get the consensus read for each drone
rule consensus:
    input:
        vcf = 'output/000_tmp/drones/{indiv}.vcf.gz',
        ref = ref
    output:
        temp('output/000_tmp/drones/{indiv}/{chr}.fa')
    log:
        'output/logs/consensus.{indiv}.{chr}.log'
    container:
        samtools
    shell:
        'samtools faidx {input.ref} {wildcards.chr} '
        '| '
        'bcftools consensus '
        '{input.vcf} '
        '>{output} '
        '2>{log}'

rule indiv_vcf:
    input:
        'output/020_filtered-genotypes/drone.renamed.vcf.gz'
    output:
        'output/000_tmp/drones/{indiv}.vcf'
    log:
        'output/logs/indiv_vcf.{indiv}.log'
    container:
        samtools
    shell:
        'bcftools view -s {wildcards.indiv} '
        '{input} '
        '> {output} '
        '2> {log}'

rule pool_vcf:
    input:
        'output/020_filtered-genotypes/pool.renamed.vcf.gz'
    output:
        temp('output/000_tmp/pools/{indiv}/{chr}.vcf')
    log:
        'output/logs/pool_vcf.{indiv}.{chr}.log'
    container:
        samtools
    shell:
        'bcftools view '
        '-s {wildcards.indiv} '
        '-r {wildcards.chr} '
        '{input} '
        '>{output} '
        '2>{log}'

# filter, split and reaheader the vcfs
rule reheader_vcf:
    input:
        vcf = 'output/020_filtered-genotypes/{set}.filtered.vcf.gz',
        headers = 'output/020_filtered-genotypes/{set}.headers.txt'
    output:
        temp('output/020_filtered-genotypes/{set}.renamed.vcf')
    log:
        'output/logs/reheader_vcf.{set}.log'
    container:
        samtools
    shell:
        'bcftools reheader '
        '-s {input.headers} '
        '{input.vcf} '
        '2>{log} '
        '| bcftools view >{output} 2>>{log} '

rule header_file:
    input:
        cnv_map
    output:
        temp('output/020_filtered-genotypes/{set}.headers.txt')
    params:
        query = lambda wildcards: f'_{wildcards.set}'
    container:
        samtools
    shell:
        'grep "{params.query}" {input} '
        '| cut -f1 '
        '| awk -F "_" \'{{print $1"_"$2" "$1}}\' '
        '> {output}'

rule split_vcf:
    input:
        vcf = 'output/020_filtered-genotypes/filtered.vcf.gz',
        cnv_map = cnv_map,
    output:
        temp('output/020_filtered-genotypes/{set}.filtered.vcf')
    log:
        'output/logs/split_vcf.{set}.log'
    params:
        query = lambda wildcards: f'_{wildcards.set}'
    container:
        samtools
    shell:
        'bcftools view '
        '-S <( grep "{params.query}" {input.cnv_map}  | cut -f1 ) '
        '{input.vcf} '
        '> {output} '
        '2> {log}'

rule filter:
    input:
        vcf = 'output/010_genotypes/calls.vcf.gz'
    output:
        temp('output/020_filtered-genotypes/filtered.vcf')
    params:
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/filter.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-m2 -M2 -v snps '  # biallelic snps only
        '--min-af {params.min_maf}:nonmajor '
        '--exclude "F_MISSING>{params.f_missing}" '
        '{input.vcf} '
        '> {output} '
        '2> {log}'

# genotype
checkpoint genotype:
    input:
        csv = sample_info,
        ref = ref,
        cnv_map = cnv_map,
        reads = expand('output/000_tmp/reads/{sample}_R{r}.fq.gz',
                       sample=all_samples,
                       r=[1, 2])
    output:
        cutoffs = 'output/010_genotypes/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/calls.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        ref = 'output/010_genotypes/015_ref/ref.fasta',
        fai = 'output/010_genotypes/015_ref/ref.fasta.fai'
    params:
        wd = 'output/010_genotypes',
    log:
        'output/logs/genotype.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--cnv_map {input.cnv_map} '
        '--threads {threads} '
        '--csd '
        '--restart_times 1 '
        '&>> {log}'


# combine reads for indivs
rule combine_reads:
    input:
        r1_1 = 'data/reads/{sample}s/{sample}s_R1.fq.gz',
        r1_2 = 'data/reads/{sample}s/{sample}s_R2.fq.gz',
        r2_1 = 'data/reads2/{sample}s/{sample}s_R1.fq.gz',
        r2_2 = 'data/reads2/{sample}s/{sample}s_R2.fq.gz'
    output:
        r1 = temp('output/000_tmp/reads/{sample}_R1.fq.gz'),
        r2 = temp('output/000_tmp/reads/{sample}_R2.fq.gz')
    singularity:
        samtools
    shell:
        'cat {input.r1_1} {input.r2_1} > {output.r1} & '
        'cat {input.r1_2} {input.r2_2} > {output.r2} & '
        'wait'



# generic bamfile subset
rule subset_bamfile:
    input:
        'output/{folder}/{file}.bam'
    output:
        'output/bam_subsets/{folder}/{file}.bam'
    log:
        'output/logs/{folder}_{file}_subset-bamfile.log'
    params:
        s = 0.001
    threads:
        2
    singularity:
        samtools
    shell:
        'samtools view '
        '-s {params.s} '
        '-O BAM '
        '{input} '
        '> {output} '
        '2> {log}'

# generic bamfile index
rule index_bamfile:
    input:
        'output/{folder}/{file}.bam'
    output:
        'output/{folder}/{file}.bam.bai'
    log:
        'output/logs/{folder}_{file}_index-bamfile.log'
    threads:
        2
    singularity:
        samtools
    shell:
        'samtools index -@ {threads} {input} 2> {log}'

# generic index rule
rule index_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    log:
        'output/logs/{folder}/{file}_index-vcf.log'
    # wildcard_constraints:
    #     file = 'filtered|phased_drones|phased_pools'
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} 2> {log} '
        '; '
        'tabix -p vcf {output.gz} 2>> {log}'


# CODE FOR DEMULTIPLEXING Undetermined FILES
# NOT ACTUALLY USEFUL
def match_failed_reads(wildcards):
    my_bc = sample_df.loc[(wildcards.run, wildcards.indiv), 'barcode']
    r1_good = sample_df.loc[(wildcards.run, wildcards.indiv), 'r1_path']
    r2_good = sample_df.loc[(wildcards.run, wildcards.indiv), 'r2_path']
    return {'r1': r1_good,
            'r2': r2_good,
            'r1_failed': f'output/000_tmp/reads/{my_bc}_r1.fastq.gz',
            'r2_failed': f'output/000_tmp/reads/{my_bc}_r2.fastq.gz'}

rule rename_target:
    input:
        expand('output/000_tmp/drones_{indiv}.fastq.gz',
               indiv=drone_samples)

rule rename_tmp:
    input:
        unpack(match_failed_reads)
    output:
        'output/000_tmp/{run}_{indiv}.fastq.gz'
    shell:
        'cat {input.r1} {input.r1_failed} > {output} ; '
        'cat {input.r2} {input.r2_failed} > {output} ; '

rule demultiplex:
    input:
        r1 = 'data/reads/failed_demultiplex/Undetermined_S0_L002_R1_001.fastq.gz',
        r2 = 'data/reads/failed_demultiplex/Undetermined_S0_L002_R2_001.fastq.gz'
    output:
        expand('output/000_tmp/reads/{barcode}_r{r}.fastq.gz',
               barcode=sorted(set(sample_df['barcode'])),
               r=['1', '2'])
    log:
        'output/logs/demultiplex.log'
    params:
        names = ','.join(sorted(set(sample_df['barcode']))),
        out = 'output/000_tmp/reads/%_r1.fastq.gz',
        out2 = 'output/000_tmp/reads/%_r2.fastq.gz'
    singularity:
        bbmap
    shell:
        'demuxbyname.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={params.out} '
        'out2={params.out2} '
        'names={params.names} '
        'prefixmode=f '
        '-Xmx100g '
        '2> {log}'


