#!/usr/bin/env nextflow

process REMOVE {
    label 'process_low'
    conda 'envs/bedtools_env.yml'
    container 'ghcr.io/bf528/bedtools:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bedA)
    path(blacklist)

    output:
    path("repr_peaks_filtered.bed"), emit: filtered

    shell:
    """
    bedtools intersect -a $bedA -b $blacklist -v > repr_peaks_filtered.bed
    """
}