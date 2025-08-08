#!/usr/bin/env nextflow

process TRIM {
    conda 'envs/trimmomatic_env.yml'
    container 'ghcr.io/bf528/trimmomatic:latest'
    publishDir params.outdir, mode: 'copy'
    
    input:
    path adapters
    tuple val(sample_id), path(reads)


    output:
    tuple val(sample_id), path("*R*P.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("*.log"), emit: log

    shell:
    """
    trimmomatic PE \\
    -trimlog ${sample_id}.log \\
    ${reads[0]} ${reads[1]} \\
    ${sample_id}_R1P.fastq.gz ${sample_id}_R1U.fastq.gz ${sample_id}_R2P.fastq.gz ${sample_id}_R2U.fastq.gz \\
    ILLUMINACLIP:${adapters}:2:30:10:2:True LEADING:3 TRAILING:3
    """

}
