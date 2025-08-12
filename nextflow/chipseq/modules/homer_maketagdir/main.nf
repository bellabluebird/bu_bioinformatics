#!/usr/bin/env nextflow

process TAGDIR {
    label 'process_medium'
    container 'ghcr.io/bf528/homer_samtools:latest'
    conda 'envs/homer_env.yml'

    input:
    tuple val(name), path(bam)

    output:
    tuple val(name), path("${name}_tags")

    script:
    """
    makeTagDirectory ${name}_tags $bam
    """
}


