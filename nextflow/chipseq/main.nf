#!/usr/bin/env nextflow

include { INTERSECT } from './modules/bedtools_intersect'
include { REMOVE } from './modules/bedtools_remove'
include { BOWTIE2_ALIGN } from './modules/bowtie2_align'
include { BOWTIE2_BUILD } from './modules/bowtie2_build'
include { BAMCOVERAGE } from './modules/deeptools_bamcoverage'
include { COMPUTEMATRIX } from './modules/deeptools_computematrix'
include { MULTIBWSUMMARY } from './modules/deeptools_multibwsummary'
include { PLOTCORRELATION } from './modules/deeptools_plotcorrelation'
include { PLOTPROFILE } from './modules/deeptools_plotprofile'
include { FASTQC } from './modules/fastqc'
include { UNZIP } from './modules/gunzip'
include { ANNOTATE } from './modules/homer_annotatepeaks'
include { FIND_MOTIFS_GENOME } from './modules/homer_findmotifsgenome'
include { FINDPEAKS } from './modules/homer_findpeaks'
include { TAGDIR } from './modules/homer_maketagdir'
include { CALLPEAKS } from './modules/macs3_callpeak'
include { MULTIQC } from './modules/multiqc'
include { RENAME_FASTQ } from './modules/rename_fastq'
include { SAMTOOLS_BAM } from './modules/samtools_bam'
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat'
include { SAMTOOLS_IDX } from './modules/samtools_idx'
include { SAMTOOLS_SORT } from './modules/samtools_sort'
include { SAMTOOLS_SORT_IDX } from './modules/samtools_sort_idx'
include { SUBSET_FASTQ } from './modules/subset_fastq'
include { SUBSET } from './modules/subset_genome'
include { TRIM } from './modules/trimmomatic'

workflow {
    Channel.fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map { row -> tuple(row.name, file(row.path)) }
        .set { sample_ch }

    // qc using fastqc
    fastqc_ch = FASTQC(sample_ch)
    
    // qc using trimmmomatic
    TRIM(sample_ch, params.adapter_fa)

    // build bowtie2 index
    BOWTIE2_BUILD(params.genome)

    // align reads to reference genome
    BOWTIE2_ALIGN(TRIM.out.trimmed, BOWTIE2_BUILD.out.index, BOWTIE2_BUILD.out.name)

    // samtools_stats - use the direct output without .sorted_bam
    SAMTOOLS_FLAGSTAT(BOWTIE2_ALIGN.out.bam)

    // combining
    combined_ch = TRIM.out.log.mix(FASTQC.out.zip).mix(SAMTOOLS_FLAGSTAT.out.flagstat).flatten().collect()

    // aggregate results with multiqc
    MULTIQC(combined_ch)
    combined_ch.view()

    // calc alignment stats using samtools flagstat
    // samtools_sort
    SAMTOOLS_SORT(BOWTIE2_ALIGN.out)
    
    // samtools_indexing - use the direct output without .sorted_bam
    SAMTOOLS_IDX(SAMTOOLS_SORT.out)

    // generate bigwig files from bam
    bam_ch = BAMCOVERAGE(SAMTOOLS_IDX.out.index)

    // bigwig matrix
    bigwig_ch = bam_ch.map{ meta, bw -> bw }.collect()
    MULTIBWSUMMARY(bigwig_ch)

    // plotting correlation
    PLOTCORRELATION(MULTIBWSUMMARY.out.multibwsummary, "pearson")

    // macs3 peaks calling; imported due to errors

    // pulling intersect from files
    Channel.fromPath("/projectnb/bf528/materials/project-2-chipseq/refs/rep*_peaks.narrowPeak")
    .collect()
    .map { files ->
            files.sort()
            tuple("repr_preaks", files[0], files[1])
    }
    .set {intersect_ch}

    // running intersect
    INTERSECT(intersect_ch)

    // Call REMOVE with the two separate inputs
    REMOVE(INTERSECT.out, params.blacklist)

    // annotating peaks to their nearest genomic feature using homer
    // For ANNOTATE, need all three inputs: peaks file, genome, and GTF
    ANNOTATE(REMOVE.out.filtered, params.genome, params.gtf)

    // running compute matrix on three inputs
    COMPUTEMATRIX(
        bam_ch,                        // BigWig files
        params.human_bed,              // Filtered peaks file
        params.window_size             // Window size parameter
    )

    // plotting a simple visualization of read counts
    PLOTPROFILE(
        COMPUTEMATRIX.out.matrix
    )

    // use findMotifsGenome.pl to perform motif enrichment analysis
    FIND_MOTIFS_GENOME(
        REMOVE.out.filtered,      // Input peaks
        params.genome             // Reference genome
    )
}
