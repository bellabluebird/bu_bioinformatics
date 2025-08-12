#!/usr/bin/env nextflow

process INTERSECT {
    label 'process_high'
    conda 'envs/bedtools_env.yml'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bedA), path(bedB)

    output:
    tuple val(meta), path("repr_peaks.bed"), emit: intersect

    shell:
    """
    bedtools intersect -a $bedA -b $bedB -f .5 -r > repr_peaks.bed
    """
}