# nextflow cheat sheet

## channels

### basic commands
```nextflow
// creation
Channel.fromPath("*.fastq")                    // individual files
Channel.fromFilePairs("*_{R1,R2}.fastq")       // paired files
Channel.of("value1", "value2")                 // static values
Channel.fromPath(params.input)                 // from config params

// transformation operators
.map{ it[x] }                    // extract x element from tuple (starts at 0)
.transpose()                     // flatten paired files: [id, [file1, file2]] → [id, file1], [id, file2]
.collect()                       // wait for all, create single list (queue → value channel)
.mix(other_channel)              // combine two channels (concatenate, not interleave)
.flatten()                       // flatten nested lists
.first()                         // take first item only
.take(5)                         // take first 5 items
.set{ my_channel }               // assign to variable

// debugging
.view()                          // print channel contents
.view{ "sample: $it" }           // custom debug message
```

### types and transformations
```nextflow
// queue channels - consumed once, multiple emissions
Channel.fromPath("*.fastq")      // ↓ file1 ↓ file2 ↓ file3

// value channels - reusable, single emission  
files.collect()                  // ↓ [file1, file2, file3] (single emission)

// common transformations
["sample1_R1.fastq", "sample1_R2.fastq"] 
→ fromFilePairs() 
→ ["sample1", [sample1_R1.fastq, sample1_R2.fastq]]

["sample1", [R1, R2]] 
→ transpose() 
→ ["sample1", R1], ["sample1", R2]

[file1, file2, file3] 
→ collect() 
→ [file1, file2, file3] (single emission for aggregation)
```

## modules and processes

### importing processes
```nextflow
include {FASTQC} from './modules/fastqc'
include {ALIGN} from './modules/star_align'
include {ALIGN as STAR_ALIGN} from './modules/star_align'    // with alias
```

### process template
```nextflow
process TOOL_NAME {
    // directives
    container 'docker/image:tag'           // or conda 'bioconda::tool'
    publishDir params.outdir, mode: "copy"
    cpus 4
    memory '8 GB'
    tag "$sample_id"                       // job identification
    when: params.run_tool                  // conditional execution
    
    input:
    tuple val(sample_id), path(reads)
    path reference
    val threshold
    
    output:
    path "${sample_id}.out", emit: result
    path "${sample_id}.log", emit: logs
    tuple val(sample_id), path("*.vcf"), emit: variants
    
    script:                               // bash (default)
    """
    tool_command -i $reads -r $reference -o ${sample_id}.out -t $threshold
    """
    
    // alternatives:
    // shell: '''command --input !{input}'''     // ! syntax for variables
    // exec: sample_id = "processed_" + sample_id // groovy code
}
```

## workflow coordination

### execution patterns
```nextflow
// parallel processing
FASTQC(sample_ch)                         // runs on all samples in parallel

// dependency handling  
INDEX(genome, gtf)                        // runs immediately
ALIGN(reads_ch, INDEX.out.index)          // waits for INDEX to complete

// collection for aggregation
PROCESS.out.files.collect() | MULTIQC     // waits for all, processes together

// channel splitting for different uses
Channel.fromFilePairs(params.reads).set { align_ch }      // paired for alignment
Channel.fromFilePairs(params.reads).transpose().set { qc_ch } // individual for qc
```

### cardinality patterns
```nextflow
// one-to-many: broadcast to all
ALIGN(reads_ch, INDEX.out.index)          // index used by all alignments

// many-to-one: collect all  
MULTIQC(all_qc_files.collect())           // waits for all qc files

// accessing process outputs
ALIGN.out.bam                             // aligned files
ALIGN.out.log                             // log files
```

### typical rnaseq flow
```nextflow
reads → FASTQC → collect → MULTIQC                    // parallel qc
reads → INDEX (reference) → ALIGN → bam files         // alignment
bam files → VERSE → count files → collect → CONCAT    // quantification
```

## execution and debugging

### running pipelines
```bash
nextflow run pipeline.nf                              // basic execution
nextflow run pipeline.nf --reads "data/*"             // override params
nextflow run pipeline.nf -resume                      // resume from failure
nextflow run pipeline.nf -profile docker              // use profile

# monitoring and reports
nextflow run pipeline.nf -with-report report.html     // execution report
nextflow run pipeline.nf -with-trace trace.txt        // process details
nextflow run pipeline.nf -with-timeline timeline.html // timeline view
```

### troubleshooting
```nextflow
// common channel error
my_ch = Channel.fromPath("*.fastq")
PROC_A(my_ch)
PROC_B(my_ch)                        // ERROR: channel consumed

// fix: create separate channels
Channel.fromPath("*.fastq").set{ ch_a }
Channel.fromPath("*.fastq").set{ ch_b }

// debugging channel flow
reads_ch.view("reads:")
ALIGN.out.bam.view("alignments:")
```

### file management
```nextflow
// work directories: temporary execution (work/xx/yyyy/)
// publishDir: final outputs copied to results/

publishDir params.outdir, mode: "copy"    // copy to final location
publishDir params.outdir, mode: "link"    // symlink to final location

// cleanup
nextflow clean -f                         // remove work directories
```

## config integration

### parameters and resources
```nextflow
// nextflow.config
params {
    reads = "data/*_{R1,R2}.fastq.gz"
    outdir = "results/"
    run_optional = false
}

process {
    withName: ALIGN {
        cpus = 8
        memory = '16 GB'
    }
}

// using in workflow
Channel.fromFilePairs(params.reads)
publishDir params.outdir
```
