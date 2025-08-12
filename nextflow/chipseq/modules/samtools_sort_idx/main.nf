#!/usr/bin/env nextflow

process SAMTOOLS_SORT_IDX {
    label 'process_single'
    conda 'envs/samtools_env.yml'
    container 'ghcr.io/bf528/samtools:latest'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.sorted.bam"), path("*.bai"), emit: index

    shell:
    """
    samtools sort -@ $task.cpus $bam > ${meta}.sorted.bam
    samtools index --threads $task.cpus ${meta}.sorted.bam
    """
}