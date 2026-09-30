process NGSPECIESID {
    tag "$barcode_id"
    label 'ngsid_env'
    publishDir { "${params.outdir}/${barcode_id}" }, mode: 'copy',
        pattern: "*.{fasta,tsv}"

    input:
    tuple val(barcode_id), path(filtered), val(medaka_model)

    output:
    tuple val(barcode_id), path("${barcode_id}_consensus_cluster*.fasta"), emit: consensus
    tuple val(barcode_id), path("${barcode_id}_cluster_read_counts.tsv"), emit: read_counts
    tuple val(barcode_id), path("${barcode_id}_cluster_sizes.tsv"),       emit: cluster_sizes

    script:
    """
    # The wrapper must win over the conda spoa: point it at the real binary, then put it first on PATH.
    export REAL_SPOA="\$(dirname "\$(command -v NGSpeciesID)")/spoa"
    export PATH="${projectDir}/bin:\$PATH"
    [ "\$(command -v spoa)" = "${projectDir}/bin/spoa" ] || { echo "spoa wrapper is not first on PATH" >&2; exit 1; }

    NGSpeciesID --ont --consensus --medaka \\
        --medaka_model ${medaka_model} \\
        --sample_size 0 \\
        --abundance_ratio ${params.abundance_ratio} \\
        --max_seqs_for_consensus ${params.max_seqs_for_consensus} \\
        --mapped_threshold ${params.mapped_threshold} \\
        --aligned_threshold ${params.aligned_threshold} \\
        --symmetric_map_align_thresholds \\
        --fastq ${filtered} \\
        --outfolder ${barcode_id}_pvcsp \\
        --t ${task.cpus}

    collect_consensus.sh ${barcode_id} ${barcode_id}_pvcsp ${filtered}
    """
}
