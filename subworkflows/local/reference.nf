//
// call tsinfer using reference alleles as ancestral alleles
//
include { TSINFER_REFERENCE } from '../../modules/local/tsinfer/reference/main'


workflow REFERENCE {
    take:
    focal_vcf_ch        // channel: focal vcf file (phased) [ meta, Path(vcf) ]
    samples_ch          // channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = channel.empty()

    // now create a tstree file
    TSINFER_REFERENCE(
        focal_vcf_ch,
        samples_ch.first()
    )
    ch_versions = ch_versions.mix( TSINFER_REFERENCE.out.versions )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
