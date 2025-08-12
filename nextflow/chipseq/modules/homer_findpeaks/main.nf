#!/usr/bin/env nextflow

process FINDPEAKS {
    label 'process_medium'
    container 'ghcr.io/bf528/homer_samtools:latest'

    input:
    tuple val(name), path(ip), path(input)

    output:
    path("${name}_peaks.txt"), emit: peaks

    script:
    """
    findPeaks $ip -style factor -o ${name}_peaks.txt -i $input
    """
}


