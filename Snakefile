#!/usr/bin/env python3

import pandas
from pathlib import Path
from tempfile import mkdtemp


def get_fastqc_reads(wildcards):
    if wildcards.type == 'raw':
        my_reads = get_reads(wildcards)
        return my_reads
    elif wildcards.type == 'processed':
        return 'output/010_process/{sample}.r1.fastq.gz'


def get_reads(wildcards):
    sample = wildcards.sample
    my_filename = sample_table.loc[sample]['filename']
    return(Path(reads_dir, my_filename).resolve().as_posix())


config_file = 'config/rnaseq_filenames.csv'
reads_dir = 'data/reads'
ref = 'data/ref/GCF_000001635.27_GRCm39_genomic.fna'
gff = 'data/ref/GCF_000001635.27_GRCm39_genomic.gff'

# containers
bbmap = ('https://github.com/deardenlab/container-bbmap/'
         'releases/download/0.0.3/container-bbmap.bbmap_38.90.sif')
fastqc = 'docker://biocontainers/fastqc:v0.11.9_cv7'
multiqc = 'docker://ewels/multiqc:1.9'
star = ('https://github.com/deardenlab/container-star/'
        'releases/download/0.0.1/container-star.2.7.9a.sif')

sample_table = pandas.read_csv(
    config_file,
    index_col="samplename")
all_samples = sorted(set(sample_table.index))


rule target:
    input:
        expand('output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
               sample=all_samples),
        'output/017_multiqc/multiqc_report.html',
        # troubleshooting GTF
        # 'output/025_star/pass1/C_M_T_61.SJ.out.tab'


rule star_second_pass:
    input:
        r1 = 'output/010_process/{sample}.r1.fastq.gz',
        star_reference = 'output/007_star-index/SA',
        junctions = expand('output/025_star/pass1/{sample}.SJ.out.tab',
                           sample=all_samples)
    output:
        bam = 'output/025_star/pass2/{sample}.Aligned.sortedByCoord.out.bam',
        counts = 'output/025_star/pass2/{sample}.ReadsPerGene.out.tab'
    threads:
        10
    params:
        genome_dir = 'output/007_star-index',
        prefix = 'output/025_star/pass2/{sample}.'
    log:
        'output/logs/star_second_pass.{sample}.log'
    resources:
        time = 59,
        mem_mb = 32 * 1000
    container:
        star
    shell:
        'STAR '
        '--runThreadN {threads} '
        '--genomeDir {params.genome_dir} '
        '--sjdbFileChrStartEnd {input.junctions} '
        '--outSAMtype BAM SortedByCoordinate '
        '--outReadsUnmapped Fastx '
        '--quantMode GeneCounts '
        '--readFilesCommand zcat '
        '--readFilesIn {input.r1} '
        '--outFileNamePrefix {params.prefix} '
        '--outTmpDir ' + mkdtemp() + ' '
        '&> {log}'


rule star_first_pass:
    input:
        r1 = 'output/010_process/{sample}.r1.fastq.gz',
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
        '--readFilesCommand zcat '
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
        '--outTmpDir ' + mkdtemp() + ' '
        '--sjdbGTFtagExonParentTranscript Parent '
        '--sjdbGTFtagExonParentGene gene '
        # '--sjdbGTFtagExonParentGeneName Name '
        '&> {log}'


# trim and qc
rule multiqc:
    input:
        expand('output/015_fastqc/{type}/{sample}.fastqc',
               type=['raw', 'processed'],
               sample=all_samples),
        expand('output/025_star/pass2/{sample}.ReadsPerGene.out.tab',
               sample=all_samples)
    output:
        'output/017_multiqc/multiqc_report.html'
    params:
        outdir = 'output/017_multiqc'
    log:
        'output/logs/multiqc.log'
    resources:
        time = 59,
    container:
        multiqc
    shell:
        'multiqc '
        '-o {params.outdir} '
        'output '
        '2> {log}'


rule fastqc:
    input:
        get_fastqc_reads
    output:
        'output/015_fastqc/{type}/{sample}.fastqc'
    params:
        outdir = 'output/015_fastqc/{type}'
    log:
        'output/logs/fastqc.{sample}.{type}.log'
    container:
        fastqc
    shell:
        'fastqc '
        '--threads {threads} '
        '-o {params.outdir} '
        '{input} '
        '&> {log} '
        '; touch {output}'


rule trim:
    input:
        get_reads
    output:
        r1 = 'output/010_process/{sample}.r1.fastq.gz'
    params:
        adapters = '/adapters.fa'
    log:
        'output/logs/trim.{sample}.log'
    threads:
        1
    resources:
        time = 25,
        mem_mb = 10 * 1000
    container:
        bbmap
    shell:
        'bbduk.sh '
        '-Xmx{resources.mem_mb}m '
        'zl=9 '
        'in={input} '
        'out={output.r1} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

