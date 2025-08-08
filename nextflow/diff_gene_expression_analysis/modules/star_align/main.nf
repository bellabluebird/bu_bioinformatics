#!/usr/bin/env nextflow

process ALIGN {
    label 'process_high'
    conda 'envs/star_env.yml'
    container 'ghcr.io/bf528/star:latest'
    publishDir params.outdir, pattern: "*.Log.final.out"

    input:
    tuple val(meta), path(reads)
    path(index)

    output:
    tuple val(meta), path("${meta}.Aligned.out.bam"), emit: bam
    tuple val(meta), path("${meta}.Log.final.out"), emit: log

    script:
    """
    STAR --runThreadN $task.cpus --genomeDir $index --readFilesIn $reads --readFilesCommand zcat --outFileNamePrefix $meta. --outSAMtype BAM Unsorted
    """
}
