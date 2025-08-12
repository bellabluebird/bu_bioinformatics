#!/usr/bin/env nextflow

process COMPUTEMATRIX {
    
    label 'process_medium'
    conda 'envs/deeptools_env.yml'
    container 'ghcr.io/bf528/deeptools:latest'

    input:
    tuple val(meta), path(bw)
    path(bed)
    val window

    output:
    tuple val(meta), path('*.gz'), emit: matrix

    shell:
    """
    computeMatrix scale-regions -S $bw -R $bed -a $window -b $window -p $task.cpus -o ${meta}_matrix.gz
    """

}