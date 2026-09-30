process DETECT_MODEL {
    tag "$barcode_id"

    input:
    tuple val(barcode_id), path(fastq)

    output:
    tuple val(barcode_id), env('MODEL')

    script:
    def override = params.medaka_model ?: ''
    """
    MODEL=\$(detect_medaka_model.sh ${fastq} '${override}')
    """
}
