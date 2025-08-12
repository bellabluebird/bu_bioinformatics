#!/usr/bin/env nextflow

process PLOTPROFILE {
    
    label 'process_low'
    conda 'envs/deeptools_env.yml'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(matrix)
    
    output:
    path('*.png'), emit: profile

    shell:
    """
    plotProfile -m $matrix -o ${meta}_signal_coverage.png
    """

}