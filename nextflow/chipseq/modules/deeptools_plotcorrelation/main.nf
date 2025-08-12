#!/usr/bin/env nextflow

process PLOTCORRELATION {
    label 'process_medium'
    conda 'envs/deeptools_env.yml'
    container 'ghcr.io/bf528/deeptools:latest'
    publishDir params.outdir

    input:
    path(multibw) 
    val(cortype)

    output:
    path("${cortype}_plot.png")

    shell:
    """
    plotCorrelation -in $multibw -c $cortype -p heatmap --plotNumbers -o ${cortype}_plot.png
    """

}






