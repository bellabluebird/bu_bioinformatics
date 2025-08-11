#!/usr/bin/env nextflow
// rnaseq pipeline for differential expression analysis!
// bella pfeiffer, boston university 2025
// see cheatsheet for more detail; reach out with questions
// i'm not an expert but hopefully this is well-annotated and easy to follow

// imports modules from individual files
include {FASTQC} from './modules/fastqc'
include {INDEX} from './modules/star_index'
include {ALIGN} from './modules/star_align'
include {SUBSET} from './modules/subset_genome'
include {VERSE} from './modules/verse'
include {CONCAT} from './modules/concat'
include {PARSE_GTF} from './modules/parse_gtf'
include {MULTIQC} from './modules/multiqc'

workflow {
    '''
    breaking channels down
        channels represent streams of data, processing in 
        parallel through the pipeline. nextflow automatically 
        distributes + coordinates our execution of processes, 
        making it much more efficient (parallelizes).

    .fromFilePairs creates tuples of file pairs
        Input files: 
            sample1_{1,2}.fastq, sample2_{1,2}.fastq
        Output: 
            ↓ ["sample1", [sample1_1.fastq, sample1_2.fastq]]
            ↓ ["sample2", [sample2_1.fastq, sample2_2.fastq]]

    .transpose() flattens paired files
        Input:  
            ["sample1", [file1, file2]]
        Output: 
            ↓ ["sample1", file1]
            ↓ ["sample1", file2]

    params.reads is a list of fastq files defined in the
    nextflow.config file, which is passed to the pipeline.
    
    .set is used to assign the channel to a variable

    fastqc_ch and align_ch are channels that will be used
    in the FASTQC and ALIGN processes respectively.

    fastqc_ch will contain the fastq files for quality control,
    while align_ch will contain the fastq files for alignment.
    The .set method allows us to reference these channels
    later in the workflow.
    '''

    Channel.fromFilePairs(params.reads).transpose().set{ fastqc_ch }
    // looks like this: 
        // ↓ ["sample1", file1]
        // ↓ ["sample1", file2]
    Channel.fromFilePairs(params.reads).set { align_ch }
    // looks like this: 
        // ↓ ["sample1", [sample1_R1.fastq, sample1_R2.fastq]]
        // ↓ ["sample2", [sample2_R1.fastq, sample2_R2.fastq]]
    // form fits function: we need paired files for alignment

    // process gtf annotation file to extract gene information
    // not sure what the output is here - need to look into it
    PARSE_GTF(params.gtf)
    
    // runs qc analysis on all fastq files
    // generates html reports + zips with metrics
    // runs all samples in parallel
    // outputs: *_fastqc.html (visual reports), *_fastqc.zip (data for multiqc)
    FASTQC(fastqc_ch)

    // creates a STAR genome index
    // combines reference genome sequence with gene annotation
    // this has nothing to do with our files, used for later analysis
    // iirc this took FOREVER to run

    // genome FASTA + GTF = splice-aware alignment
    // output: INDEX.out.INDEX - directory containing STAR index files
    INDEX(params.genome, params.gtf)

    // takes paired fastq files + genome index from INDEX
    // aligns rna-seq reads to reference genome + handles splice junctions
    // produces aligned reads in BAM
    // outputs: ALIGN.out.bam (aligned reads), ALIGN.out.log (alignment stats for multiqc)
    ALIGN(align_ch, INDEX.out.index)

    // ** NEXT SECTION ** 
    // collecting and preparing outputs for multiqc

    // .FASTQC.out.zip: gets the fastqc output files
        // ↓ ["sample1_R1", "/path/sample1_R1_fastqc.zip"]
    // .map.it{ it[1] }: maps the output to get the second element of each tuple
        // ↓ "/path/sample1_R1_fastqc.zip"
    // .collect(): waits for all fastqc outputs to be ready,
    // and collects them into a value channel called fastqc_out
        // ↓ ["/path/sample1_R1_fastqc.zip", 
        //   ... "/path/sample3_R2_fastqc.zip"]
    FASTQC.out.zip.map{ it[1] }.collect()
    | set { fastqc_out }

    // .ALIGN.out.log: gets the alignment log files
    // .map.it{ it[1] }: maps the output to get the second element
    // of each tuple, which is the log file path
    // .collect(): waits for all alignment logs to be ready,
    // and collects them into a value channel called star_log
    ALIGN.out.log.map{ it[1] }.collect()
    | set { star_log }

    // this one is a little more complex
    // starting channels: 
        // fastqc_out: ↓ [fastqc_zip1, fastqc_zip2, fastqc_zip3, 
        //                fastqc_zip4, fastqc_zip5, fastqc_zip6]
        // star_log:   ↓ [star_log1, star_log2, star_log3]
    // .mix(star_log): combines the fastqc outputs with the star logs
    // .collect(): waits for all mixed outputs to be ready,
    // and collects them into a value channel called multiqc_ch
    // output:
        // ↓ [fastqc_zip1, fastqc_zip2, fastqc_zip3, fastqc_zip4, 
        //    fastqc_zip5, fastqc_zip6, star_log1, star_log2, star_log3]
    fastqc_out.mix(star_log).collect()
    | set { multiqc_ch }

    // input: single list with all qc files
    // output: unified html report for multiqc
    // the report will be saved in the current working directory
    // BECAUSE of the above collect statements,
    // we know multiqc_ch will contain all the files
    MULTIQC(multiqc_ch)

    // input: ALIGN.out.bam (aligned reads in BAM format)
        // ↓ ["sample1", "/path/sample1_Aligned.sortedByCoord.out.bam"]
    // input: params.gtf (gene annotation file)
    // output: VERSE.out.counts (gene expression counts)
        // ↓ ["sample1", "/path/sample1.counts.txt"]
        //   TXT: ENSG00000000003    1205

    // VERSE is a tool for quantifying gene expression
    // it reads aligned BAM files and counts reads per gene using
    // the provided GTF annotation file to define gene boundaries
    VERSE(ALIGN.out.bam, params.gtf)

    // input: VERSE.out.counts (gene expression counts)
    // output: concat_ch (channel containing all count files in one list)
    VERSE.out.counts.map{ it[1] }.collect().set{ concat_ch } 

    // input: concat_ch (channel with all count files)
    // output: count matrix in TSV format
    CONCAT(concat_ch)

    // show html report in the browser
    multiqc_ch.view()
}
