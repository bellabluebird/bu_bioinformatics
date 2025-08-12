#!/usr/bin/env nextflow

process SAMTOOLS_BAM {
    label 'process_low'
    conda 'envs/samtools_env.yml'
    container 'ghcr.io/bf528/samtools:latest'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.bam")

    shell:
    """
    samtools view -bS $sam > ${meta}.bam
    """
}