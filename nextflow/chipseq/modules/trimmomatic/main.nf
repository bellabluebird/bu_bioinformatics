#!/usr/bin/env nextflow

process TRIM {
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir params.outdir, mode: 'copy'
    
    input:
    tuple val(sample_id), path(fastq)
    path(adapters)
    
    output:
    path("*.log"), emit: log
    tuple val(sample_id), path("${fastq.simpleName}_trimmed.fastq.gz"), emit: trimmed

    shell:
    """
    trimmomatic SE $fastq ${fastq.simpleName}_trimmed.fastq.gz ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 2> ${fastq.simpleName}_trim.log
    """

}
