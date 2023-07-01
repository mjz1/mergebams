#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2

params.csv = null
params.samples = null
params.bams = null
params.barcodes = null
params.labels = null
params.outdir = "./output"


log.info """\
    M E R G E B A M S   P I P E L I N E
    ===================================
    input_csv    : ${params.csv}
    samples      : ${params.samples}
    inbams       : ${params.bams}
    barcodes     : ${params.barcodes}
    labels       : ${params.labels}
    outdir       : ${params.outdir}
    """
    .stripIndent()



process MERGEBAMS {
    label 'multicore'
    publishDir path: params.outdir, mode:'copyNoFollow'

    tag "Sample: $sampleid"

    input: 
        // tuple val(sampleid), path(bams) val(labels)
        tuple val(sampleid), path(bams, stageAs: "?/*"), val(labels), path(barcodes, stageAs: "?/*")
    
    output:
        path("${sampleid}/*")

    script:
    """
    echo "SampleID: ${sampleid}"
    echo "BAMs: ${bams}"
    echo "Labels: ${labels}"
    echo "Barcodes: ${barcodes}"

    outdir=${sampleid}
    mkdir -p \$outdir

    # Create comma sep strings for the three variables
    bam_in=\$(echo "${bams}" | tr ' ' ',')
    barcodes_in=\$(echo "${barcodes}" | tr ' ' ',')

    mergebams -i \${bam_in} -l ${labels} -b \${barcodes_in} -o \$outdir -t $task.cpus
    """
}

workflow {
    
    if (params.csv) {
        input_ch = Channel.fromPath(params.csv).
            splitCsv(header:true).
            map(row -> tuple "${row.sample_id}", "${row.bam}", "${row.label}", "${row.barcodes}")
    } else {
        input_ch = Channel.fromList(params.samples).
            merge(Channel.fromList(params.bams)).
            merge(Channel.fromList(params.labels)).
            merge(Channel.fromList(params.barcodes))
    }

    input_ch = input_ch.
            groupTuple(by: [0]).
            map(row -> tuple row[0], row[1], row[2].join(","), row[3])

    MERGEBAMS(input_ch)
}

