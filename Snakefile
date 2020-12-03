#!/usr/bin/env python3

import multiprocessing
import pandas
import tempfile


###########
# GLOBALS #
###########

bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.6')
minimap = 'shub://TomHarrop/align-utils:minimap2_2.17r941'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
whatshap = 'shub://TomHarrop/variant-utils:whatshap_491ec8e'
phase_honeybee_vcf = 'shub://TomHarrop/phase-honeybee-vcf:phase_honeybee_vcf_v0.0.2'

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
# both_indivs = ['BB34', 'BB42'] # for testing

# read reference data
fai_pd = pandas.read_csv(fai, sep='\t', header=None)
all_chr = sorted(set(fai_pd[0]))
autosomes = [x for x in all_chr if x.startswith('NC_')]


#########
# RULES #
#########

wildcard_constraints:
    chr = '|'.join(autosomes)

# rules
rule target:
    input:
        'output/045_phased/autosomes.vcf.gz'

rule aio_phase:
    input:
        ref = ref,
        vcf = 'output/020_filtered-genotypes/filtered.vcf.gz',
        bam = 'output/010_genotypes/merged.bam',
        csv = 'data/phasing.csv'
    output:
        'output/047_phasing/phased.vcf.gz'
    log:
        'output/logs/aio_phase.log'
    params:
        outdir = 'output/047_phasing'
    threads:
        workflow.cores
    container:
        phase_honeybee_vcf
    shell:
        'phase_honeybee_vcf '
        '--threads {threads} '
        '--ref {input.ref} '
        '--vcf {input.vcf} '
        '--bam {input.bam} '
        '--samples_csv {input.csv} '
        '--outdir {params.outdir} '
        '2> {log}'


rule concat:
    input:
        expand('output/040_phased-chrs/{chr}.vcf.gz',
               chr=autosomes)
    output:
        'output/045_phased/autosomes.vcf'
    log:
        'output/logs/concat.log'
    container:
        samtools
    shell:
        'bcftools concat '
        '{input} '
        '>>{output} '
        '2> {log}'

rule phase:
    input:
        pool_vcf = 'output/035_renamed_vcfs/pools/{chr}.vcf.gz',
        pool_bam = 'output/037_merged-bams/pools/{chr}.bam',
        pool_bai = 'output/037_merged-bams/pools/{chr}.bam.bai',
        drone_bam = 'output/037_merged-bams/drones/{chr}.bam',
        drone_bai = 'output/037_merged-bams/drones/{chr}.bam.bai',
        ref = ref
    output:
        temp('output/040_phased-chrs/{chr}.vcf')
    log:
        'output/logs/phase.{chr}.log'
    container:
        whatshap
    shell:
        'whatshap phase '
        '--reference {input.ref} '
        '-o {output} '
        '{input.pool_vcf} '
        '{input.drone_bam} '
        '{input.pool_bam} '
        '&>{log}'

# merge all the indivs into a per-chromsome bam / vcf
rule merge_bam:
    input:
        expand('output/000_tmp/{{type}}/{indiv}/{{chr}}.bam',
               indiv=both_indivs)
    output:
        'output/037_merged-bams/{type}/{chr}.bam'
    log:
        'output/logs/merge_bam.{type}.{chr}.log'
    container:
        samtools
    shell:
        'samtools merge {output} {input} 2>{log}'


rule sort_vcf:
    input:
        'output/000_tmp/pools/{chr}.vcf'
    output:
        temp('output/035_renamed_vcfs/pools/{chr}.vcf')
    log:
        'output/logs/sort_vcf.{chr}.log'
    container:
        samtools
    shell:
        'bcftools sort '
        '--temp-dir ' + tempfile.mkdtemp() + ' '
        '{input} '
        '>{output} '
        '2>{log}'

rule merge_vcf:
    input:
        expand('output/000_tmp/pools/{indiv}/{{chr}}.vcf.gz',
               indiv=both_indivs)
    output:
        temp('output/000_tmp/pools/{chr}.vcf')
    log:
        'output/logs/merge_vcf.{chr}.log'
    container:
        samtools
    shell:
        'bcftools merge '
        '-O v '
        '{input} '
        '>>{output} '
        '2>{log}'

# rename the pool reads
rule rename_bam:
    input:
        'output/000_tmp/pools/{indiv}/{chr}.sam'
    output:
        temp('output/000_tmp/pools/{indiv}/{chr}.bam')
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
        bam = temp('output/000_tmp/drones/{indiv}/{chr}.bam'),
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
        temp('output/000_tmp/drones/{indiv}/{chr}.renamed.fa')
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
        temp('output/000_tmp/drones/{indiv}.vcf')
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
