#!/usr/bin/env nextflow

process VERSE {
    label 'process_single'
    container 'ghcr.io/bf528/verse:latest'
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(meta), path(bam)
    path(gtf)

    output:
    tuple val(meta), path("*exon.txt"), emit: counts

    shell:
    """
    verse -S -a $gtf -o $meta $bam
    """
}
