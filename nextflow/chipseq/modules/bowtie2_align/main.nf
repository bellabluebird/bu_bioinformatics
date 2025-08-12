#!/usr/bin/env nextflow

process BOWTIE2_ALIGN {
    label 'process_high'
    container 'ghcr.io/bf528/bowtie2:latest'
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(meta), path(fastq)
    path bt2
    val name

    output:
    tuple val(meta), path('*.bam'), emit: bam

    shell:
    """
    bowtie2 -p 8 -x $bt2/$name -U $fastq | samtools view -bS > ${meta}.bam
    """
}