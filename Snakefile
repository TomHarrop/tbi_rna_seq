#!/usr/bin/env python3

import pandas
from pathlib import Path


def get_reads(wildcards):
    sample = wildcards.sample
    my_filename = sample_table.loc[sample]['filename']
    return(Path(reads_dir, my_filename).resolve().as_posix())


config_file = 'config/sample_table.csv'
reads_dir = 'data/reads'
ref = 'data/ref/GCF_000001635.27_GRCm39_genomic.fna'
gff = 'data/ref/GCF_000001635.27_GRCm39_genomic.gff'

# containers
bbmap = 'docker://ghcr.io/deardenlab/container-bbmap:bbmap_38.90'
star = 'shub://TomHarrop/align-utils:star_2.7.6a'

sample_table = pandas.read_csv(
    config_file,
    index_col="sample_name")
all_samples = sorted(set(sample_table.index))


rule target:
    input:
        expand('output/025_star/pass1/{sample}.SJ.out.tab',
               sample=all_samples)


rule star_first_pass:
    input:
        r1 = 'output/010_process/{sample}.r1.fastq',
        star_reference = 'output/007_star-index/SA'
    output:
        sjdb = 'output/025_star/pass1/{sample}.SJ.out.tab'
    threads:
        10
    params:
        genome_dir = 'output/007_star-index',
        prefix = 'output/025_star/pass1/{sample}.'
    log:
        'output/logs/star_first_pass.{sample}.log'
    resources:
        time = 59,
        mem_mb = 32 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--outSJfilterReads Unique '
        '--outSAMtype None '          # troubleshoot gtf
        # '--outSAMtype SAM '               # troubleshoot gtf
        # '--quantMode GeneCounts '       # troubleshoot gtf
        '--readFilesIn {input.r1} '
        '--outFileNamePrefix {params.prefix} '
        '&> {log}'


rule star_index:
    input:
        fasta = ref,
        gff = gff
    output:
        'output/007_star-index/SA'
    params:
        outdir = 'output/007_star-index'
    log:
        'output/logs/star_index.log'
    threads:
        10
    resources:
        time = 30,
        mem_mb = 32 * 1000
    container:
        star
    shell:
        'STAR '
        'runThreadN {threads} '
        '--runMode genomeGenerate '
        '--genomeDir {params.outdir} '
        '--genomeFastaFiles {input.fasta} '
        '--sjdbGTFfile {input.gff} '
        '--genomeSAindexNbases 12 '
        '--sjdbGTFtagExonParentTranscript Parent '
        '--sjdbGTFtagExonParentGene locus_tag '
        # '--sjdbGTFtagExonParentGeneName Name '
        '&> {log}'


rule trim:
    input:
        get_reads
    output:
        r1 = 'output/010_process/{sample}.r1.fastq'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/trim.{sample}.log'
    threads:
        1
    resources:
        time = 20,
    container:
        bbmap
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out={output.r1} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

