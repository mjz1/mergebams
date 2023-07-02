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
    publishDir path: params.outdir, mode:'copyNoFollow', overwrite: true
    cache 'lenient'
    cpus 1
    memory '32 GB'

    tag "Sample: $sampleid"

    input: 
        // tuple val(sampleid), path(bams) val(labels)
        tuple val(sampleid), path(bams, stageAs: "?/*"), val(labels), path(barcodes, stageAs: "?/*")
    
    output:
        tuple val(sampleid), path("${sampleid}/out_bam.bam"), path("${sampleid}/fail_bam.bam"), path("${sampleid}/out_barcodes.tsv.gz")

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

process SAMTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_medium'
    publishDir path: params.outdir, mode:'copyNoFollow', overwrite: true

    conda "bioconda::samtools=1.17"
    container "quay.io/biocontainers/samtools:1.17--hd87286a_1"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}/*.bam*"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    outdir=${prefix}
    mkdir -p \$outdir
    samtools sort \\
        $args \\
        -@ $task.cpus \\
        -o \$outdir/${prefix}.bam \\
        -T $prefix \\
        $bam

    samtools index \\
        $args \\
        -@ $task.cpus \\
        \$outdir/${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
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

    sam_inch = MERGEBAMS(input_ch).map {
        sample_id, bam, fbam, bc ->
            tuple([id: sample_id, single_end:false], bam)
    }

    SAMTOOLS_SORT_INDEX(sam_inch)

}

