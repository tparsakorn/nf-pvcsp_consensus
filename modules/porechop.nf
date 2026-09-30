process PORECHOP {
    tag "$barcode_id"
    label 'ngs_env'

    input:
    tuple val(barcode_id), path(fastq)

    output:
    tuple val(barcode_id), path("${barcode_id}.trim.fastq"), emit: trimmed

    script:
    """
    porechop_abi -i ${fastq} -o ${barcode_id}.trim.fastq --threads ${task.cpus}
    """
}
