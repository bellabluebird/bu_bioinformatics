#!/usr/bin/env nextflow

process ANNOTATE {
    label 'process_low'
    conda 'envs/homer_env.yml'
    container 'ghcr.io/bf528/homer:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    path(repr_filtered_peaks)
    path(genome)
    path(gtf)

    output:
    path('annotated_peaks.txt')

    shell:
    """
    annotatePeaks.pl $repr_filtered_peaks $genome -gtf $gtf > annotated_peaks.txt
    """
}



