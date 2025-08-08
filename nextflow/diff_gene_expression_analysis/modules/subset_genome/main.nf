#!/usr/bin/env nextflow

process SUBSET {
    conda 'envs/biopython_env.yml'
    container 'ghcr.io/bf528/biopython:latest'
    publishDir params.refdir, mode: 'copy'
    
    input:
    path genome
    val chrom

    output:
    path("${chrom}_only.fa"), emit: subset

    script:
    """
    subset_chrom.py -i $genome -o ${chrom}_only.fa -s ${chrom}
    """

}
