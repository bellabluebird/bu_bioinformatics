#!/usr/bin/env nextflow

process SUBSET_FASTQ {
    conda 'envs/seqtk_env.yml'
    publishDir params.refdir, mode: 'copy'
    
    input:
    path(fastq)
    val(subsample)

    output:
    path("*.subset.fastq.gz")

    script:
    """
    seqtk sample -s 1337 $fastq $subsample > ${fastq.simpleName}.subset.fastq
    gzip ${fastq.simpleName}.subset.fastq
    """
}