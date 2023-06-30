#!/usr/bin/env nextflow

// Declare syntax version
nextflow.enable.dsl=2


log.info """\
    M E R G E B A M S   P I P E L I N E
    ===================================
    samples      : ${params.samples}
    inbams       : ${params.bams}
    barcodes     : ${params.barcodes}
    labels       : ${params.labels}
    outdir       : ${params.outdir}
    """
    .stripIndent()



process MERGEBAMS {
    publishDir path: params.outdir, mode:'copy'

    input: 
        // tuple val(sampleid), path(bams) val(labels)
        tuple val(sampleid), path(bams), val(labels), path(barcodes)
    
    output:
        path(sampleid)

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

    mergebams -i \${bam_in} -l ${labels} -b \${barcodes_in} -o \$outdir
    """
}

workflow {
    input_ch = Channel.fromList(params.samples).
        merge(Channel.fromList(params.bams)).
        merge(Channel.fromList(params.labels).map{element -> element + "_"}).
        merge(Channel.fromList(params.barcodes)).
        groupTuple(by: [0]).
        map(row -> tuple row[0], row[1], row[2].join(","), row[3])

    MERGEBAMS(input_ch)
}

