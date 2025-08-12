#!/usr/bin/env nextflow

process FASTQC {
    label 'process_low'
    conda 'envs/fastqc_env.yml'
    container 'ghcr.io/bf528/fastqc:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html")
    path("*.zip"), emit: zip

    shell:
    """
    fastqc -t $task.cpus ${reads}
    """
}