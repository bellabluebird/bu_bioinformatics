#!/usr/bin/env nextflow

process MULTIBWSUMMARY {
    label 'process_medium'
    conda 'envs/deeptools_env.yml'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir

    input:
    path(bws)

    output:
    path('bw_all.npz'), emit: multibwsummary

    shell:
    """
    multiBigwigSummary bins -b  ${bws.join(' ')} --labels ${bws.baseName.join(' ')} -o bw_all.npz -p $task.cpus
    """

}