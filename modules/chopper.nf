process CHOPPER {
    tag "$barcode_id"
    label 'ngs_env'

    input:
    tuple val(barcode_id), path(trimmed)

    output:
    tuple val(barcode_id), path("${barcode_id}.filter.fastq"), emit: filtered

    script:
    """
    chopper -q ${params.min_quality} -l ${params.min_length} --maxlength ${params.max_length} \\
        --threads ${task.cpus} -i ${trimmed} > ${barcode_id}.filter.fastq
    """
}
