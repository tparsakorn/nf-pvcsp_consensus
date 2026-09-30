include { NANOPLOT     } from './modules/nanoplot'
include { PORECHOP     } from './modules/porechop'
include { CHOPPER      } from './modules/chopper'
include { DETECT_MODEL } from './modules/detect_model'
include { NGSPECIESID  } from './modules/ngspeciesid'
include { CHECK_ORF    } from './modules/check_orf'

workflow {
    if (!params.samplesheet) {
        error "Missing --samplesheet. Provide a CSV with columns barcode_id,fastq (see assets/samplesheet_example.csv)"
    }

    samples = channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> tuple(row.barcode_id, file(row.fastq, checkIfExists: true)) }

    NANOPLOT(samples)

    PORECHOP(samples)
    CHOPPER(PORECHOP.out.trimmed)

    // medaka model is read from the RAW fastq header (RG tag), unless --medaka_model is given
    DETECT_MODEL(samples)

    NGSPECIESID(CHOPPER.out.filtered.join(DETECT_MODEL.out))

    CHECK_ORF(NGSPECIESID.out.consensus)
}
