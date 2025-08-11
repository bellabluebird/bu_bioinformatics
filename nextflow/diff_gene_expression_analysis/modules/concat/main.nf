#!/usr/bin/env nextflow

process CONCAT {
    '''
    CONCAT process for merging individual gene count files into a count matrix
    
    container uses pandas image since we need python + pandas for data manipulation
    the concat_df.py script requires pandas to read and merge the count files
    
    publishDir copies the final count matrix to the results directory
    this is the end product we want to save - the analysis-ready matrix
    
    input takes a list of count file paths from the collected channel
    this comes from: VERSE.out.counts.map{ it[1] }.collect()
    which gives us all the individual count files in one list
    
    output creates the merged count matrix as a csv file
    emit: concat allows us to reference this as CONCAT.out.concat later
    
    script calls the python script we analyzed earlier
    -i flag passes the list of input files
    -o flag specifies the output filename
    the script reads each count file, extracts sample names from filenames,
    and merges them into a genes x samples matrix
    '''
    container 'ghcr.io/bf528/pandas:latest'
    publishDir params.outdir, mode: "copy"
    // copies final count matrix to results directory

    input:
    path df_list
    // collected list of count files:
        // ↓ [sample1.counts.txt, sample2.counts.txt, sample3.counts.txt]

    output:
    path 'verse_concat.csv', emit: concat
    // final count matrix:
        // ↓ verse_concat.csv (genes as rows, samples as columns)

    script:
    """
    concat_df.py -i $df_list -o verse_concat.csv
    """
    // runs python script to merge individual count files
    // creates analysis-ready count matrix for downstream tools like DESeq2
}
