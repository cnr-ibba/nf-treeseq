//
// call tsinfer using reference alleles as ancestral alleles
//
include { TSINFER_REFERENCE } from '../../modules/local/tsinfer_reference'


workflow REFERENCE {
    take:
    focal_vcf_ch        // Channel: focal vcf file (phased) [ meta, Path(vcf) ]
    samples_ch          // Channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = Channel.empty()

    // now create a tstree file
    TSINFER_REFERENCE(
        focal_vcf_ch,
        samples_ch.first()
    )

    emit:

    versions       = ch_versions                    // channel: [ versions.yml ]
}
