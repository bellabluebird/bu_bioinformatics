#!/usr/bin/env nextflow

process CALLPEAKS {
    label 'process_medium'
    conda 'envs/macs3_env.yml'
    container 'ghcr.io/bf528/macs3:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(rep), path(IP), path(CONTROL)
    val(macs3_genome)

    output:
    tuple val(rep), path('*.narrowPeak'), emit: peaks

    shell:
    """
    macs3 callpeak -t $IP -c $CONTROL -f BAM -g $macs3_genome -n $rep
    """
}