#!/usr/bin/env nextflow

process CONCAT {
    container 'ghcr.io/bf528/pandas:latest'
    publishDir params.outdir, mode: "copy"

    input:
    path df_list

    output:
    path 'verse_concat.csv', emit: concat

    script:
    """
    concat_df.py -i $df_list -o verse_concat.csv
    """

}
