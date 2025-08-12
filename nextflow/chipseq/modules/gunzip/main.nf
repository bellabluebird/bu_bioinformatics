#!/usr/bin/env nextflow

process UNZIP {
    label 'process_single'

    input:
    path gz
    
    output:
    path gz.baseName, emit: unzip

    shell:
    """
    gunzip -c ${gz} > ${gz.baseName}
    """

}