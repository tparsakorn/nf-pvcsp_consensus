process NANOPLOT {
    tag "$barcode_id"
    label 'ngs_env'
    publishDir { "${params.outdir}/${barcode_id}" }, mode: 'copy'

    input:
    tuple val(barcode_id), path(fastq)

    output:
    tuple val(barcode_id), path("nanoplot")

    script:
    // --no_static: HTML report + stats only. Static PNGs need kaleido, which needs a Chrome install on some setups.
    """
    NanoPlot --fastq ${fastq} --outdir nanoplot --threads ${task.cpus} --no_static
    """
}
