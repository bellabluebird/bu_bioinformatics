#!/usr/bin/env nextflow

process RENAME_FASTQ {
    container 'ghcr.io/bf528/biopython:latest'
    publishDir params.refdir, mode: 'move'
    
    input:
    tuple path(fastq), val(filename), val(samplename), val(new_name)

    output:
    path("*.fastq.gz")

    script:
    """
    subset_fastq.py -i $fastq -o ${new_name}.fastq -s $samplename -n $new_name
    gzip ${new_name}.fastq
    """

}