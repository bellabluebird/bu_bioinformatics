#!/usr/bin/env nextflow

process SAMTOOLS_SORT {
    label 'process_single'
    conda 'envs/samtools_env.yml'
    container 'ghcr.io/bf528/samtools:latest'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*sorted.bam")

    shell:
    """
    samtools sort -@ $task.cpus $bam > ${meta}.sorted.bam
    """
}