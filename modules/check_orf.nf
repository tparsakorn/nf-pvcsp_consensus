process CHECK_ORF {
    tag "$barcode_id"
    label 'ngs_env'
    publishDir { "${params.outdir}/${barcode_id}" }, mode: 'copy'

    input:
    tuple val(barcode_id), path(consensus)

    output:
    tuple val(barcode_id), path("${barcode_id}_check_orf.tsv")

    script:
    """
    # exit code 1 just means "at least one FAIL"; the verdict is in the report, so don't abort the run
    check_orf.py ${consensus} > ${barcode_id}_check_orf.tsv || true
    """
}
