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
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.3'
    '@e7e37748bde42ab8d6ad8dffecd5ca008089276c')
r = ('shub://TomHarrop/r-containers:r_3.6.1'
     '@e1eb426cd153fd0669bc24508673228d2f25dd76')
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'       # fixme
vcftools = ('shub://TomHarrop/variant-utils:vcftools_0.1.16'
            '@d64cc5a37951760be575c43024c66e69b2563166')
whatshap = 'shub://TomHarrop/variant-utils:whatshap_0.18'

ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
drones_csv = 'data/drones_indivonly.csv'
pools_csv = 'data/pools_indivonly.csv'

########
# MAIN #
########

# map barcodes to names
pool_data = pandas.read_csv(pools_csv,
                            index_col='sample')
pool_samples = sorted(set(pool_data.index))
drone_data = pandas.read_csv(drones_csv,
                             index_col='sample')
drone_samples = sorted(set(drone_data.index))

# generate output files for demultiplexing
sample_data = {'drones': drone_data,
               'pools': pool_data}

sample_df = pandas.concat(sample_data)

#########
# RULES #
#########

# rules
rule target:
    input:
        expand('output/020_filtered-genotypes/{run}/filtered.vcf.gz',
               run=['pools', 'drones']),
        expand('output/040_phased-indivs/{indiv}.gtf',
               indiv=pool_samples),
        'output/030_phased/phased_pools.vcf.gz',
        'output/030_phased/phased_pools.bam.bai'

rule phase_reads:
    input:
        vcf = 'output/030_phased/phased_pools.vcf.gz',
        ref = 'output/010_genotypes/pools/015_ref/ref.fasta',
        bam = 'output/010_genotypes/pools/merged.bam'
    output:
        bam = 'output/030_phased/phased_pools.bam'
    log:
        'output/logs/phase_reads.log'
    singularity:
        whatshap
    shell:
        'whatshap haplotag '
        '-o {output.bam} '
        '--reference {input.ref} '
        '{input.vcf} '
        '{input.bam} '
        '&> {log} '
        '|| true '  # don't know why it's failing but i need to see the output

rule gtf_blocks:
    input:
        vcf = 'output/030_phased/phased_pools.vcf.gz',
        fai = 'output/010_genotypes/pools/015_ref/ref.fasta.fai'
    output:
        gtf = 'output/040_phased-indivs/{indiv}.gtf',
        stats = 'output/040_phased-indivs/{indiv}.txt',
    log:
        'output/logs/040_phased-indivs/{indiv}_gtf-blocks.log'
    singularity:
        whatshap
    shell:
        'whatshap stats '
        '--gtf {output.gtf} '
        '--chr-lengths <(cut -f1,2 {input.fai}) '
        '--sample {wildcards.indiv} '
        '{input.vcf} '
        '> {output.stats} '
        '2> {log}'


rule phase_pools:
    input:
        vcf = 'output/020_filtered-genotypes/pools/filtered.vcf.gz',
        bam = 'output/010_genotypes/pools/merged.bam',
        ref = 'output/010_genotypes/pools/015_ref/ref.fasta',
        fai = 'output/010_genotypes/pools/015_ref/ref.fasta.fai',
        drone_vcf = 'output/030_phased/phased_drones.vcf.gz'
    output:
        'output/030_phased/phased_pools.vcf'
    log:
        'output/logs/phase_pools.log'
    singularity:
        whatshap
    shell:
        'whatshap phase '
        '--reference {input.ref} '
        '-o {output} '
        '--indels '
        '{input.vcf} '
        '{input.bam} '
        '{input.drone_vcf} ' # add phased drone vcf here
        '&> {log}'

# manually add a dummy phase set block to all drone calls
rule add_phase_set:
    input:
        vcf = 'output/020_filtered-genotypes/drones/filtered.vcf'
    output:
        vcf = 'output/030_phased/phased_drones.vcf',
        tmp = temp('output/030_phased/phased_drones_noheader.vcf')
    log:
        'output/logs/add_phase_set.log'
    singularity:
        r
    shell:
        'Rscript src/manual_vcf_parse.R {input.vcf} {output.tmp} '
        '&> {log} ; '
        'cat <(grep "^##" {input.vcf}) '
        '<(echo \'##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\') '
        '<(grep "^#[^#]" {input.vcf}) '
        ' {output.tmp} '
        ' > {output.vcf} '

rule filter:
    input:
        cutoffs = 'output/010_genotypes/{run}/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/{run}/calls.vcf.gz'
    output:
        'output/020_filtered-genotypes/{run}/filtered.vcf'
    params:
        min_depth = get_min_cutoff,
        max_depth = get_max_cutoff,
        maf = 0.1,
        max_missing = 0.9,
        qual = 30
    log:
        'output/logs/filter_{run}.log'
    singularity:
        vcftools
    shell:
        'vcftools '
        '--gzvcf {input.vcf} '
        '--maf {params.maf} '
        '--max-missing {params.max_missing} '
        '--minQ {params.qual} '
        '--min-meanDP {params.min_depth} '
        '--max-meanDP {params.max_depth} '
        '--minDP {params.min_depth} '
        '--maxDP {params.max_depth} '
        '--recode '
        '--stdout '
        '> {output} '
        '2> {log}'

# genotype
checkpoint genotype:
    input:
        csv = 'data/{run}_indivonly.csv',
        ref = ref
    output:
        cutoffs = 'output/010_genotypes/{run}/040_stats/ldepth.mean_cutoffs.csv',
        vcf = 'output/010_genotypes/{run}/calls.vcf.gz',
        bam = 'output/010_genotypes/{run}/merged.bam',
        ref = 'output/010_genotypes/{run}/015_ref/ref.fasta',
        fai = 'output/010_genotypes/{run}/015_ref/ref.fasta.fai'
    params:
        wd = 'output/010_genotypes/{run}',
        ploidy = lambda wildcards:
            '1' if wildcards.run == 'drones' else '2'
    log:
        'output/logs/genotype_{run}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--ploidy {params.ploidy} '
        '--threads {threads} '
        '--restart_times 1 '
        '&>> {log}'

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
    wildcard_constraints:
        file = 'filtered|phased_drones|phased_pools'
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
