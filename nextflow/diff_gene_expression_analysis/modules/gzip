#!/usr/bin/env nextflow

process ZIP {
    label 'process_single'
    conda 'envs/seqtk_env.yml'
    publishDir params.refdir

    input:
    path(fastq)
    
    output:
    path("*.gz")

    shell:
    """
    gzip $fastq
    """

}
