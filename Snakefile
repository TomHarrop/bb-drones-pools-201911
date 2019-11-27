#!/usr/bin/env python3

import multiprocessing
from pathlib import Path


def resolve_path(x):
    return str(Path(x).resolve())


honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.1'
    '@05814804f129b65b889723bf4bd7c9d40c68fe1a')
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'       # fixme
vcftools = ('shub://TomHarrop/variant-utils:vcftools_0.1.16'
            '@d64cc5a37951760be575c43024c66e69b2563166')


ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
drones_csv = 'data/drones.csv'
pools_csv = 'data/pools.csv'

# vcftools nasty output, have to run once each


ext_to_arg = {
    'frq': 'freq2 --max-alleles 2',
    'idepth': 'depth',
    'ldepth.mean': 'site-mean-depth',
    'lqual': 'site-quality',
    'imiss': 'missing-indv',
    'lmiss': 'missing-site'}


# rules
rule target:
    input:
        expand('output/010_genotypes/{run}/calls.vcf.gz',
               run=['pools', 'drones']),
        expand('output/020_stats/{run}/stats.{ext}',
               ext=list(ext_to_arg.keys()),
               run=['pools', 'drones'])

# qc
rule vcf_stats:
    input:
        'output/010_genotypes/{run}/calls.vcf.gz'
    output:
        'output/020_stats/{run}/stats.{ext}'
    log:
        'output/logs/020_stats/{run}-{ext}.log'
    params:
        wd = 'output/020_stats/{run}',
        arg = lambda wildcards: ext_to_arg[wildcards.ext]
    singularity:
        vcftools
    shell:
        'cd {params.wd} || exit 1 ; '
        'vcftools '
        '--gzvcf '
        + resolve_path('{input}') + ' '
        '--{params.arg} '
        '--out stats '
        '2> ' + resolve_path('{log}')

# genotype
rule genotype:
    input:
        csv = 'data/{run}.csv',
        ref = ref
    output:
        'output/010_genotypes/{run}/calls.vcf.gz',
        'output/010_genotypes/{run}/merged.bam'
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
        '&> {log}'


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
