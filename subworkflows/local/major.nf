//
// call tsinfer using major alleles as ancestral alleles
//
include { TSINFER_MAJOR } from '../../modules/local/tsinfer_major'


workflow MAJOR {
    take:
    focal_vcf_ch        // Channel: focal vcf file (phased) [ meta, Path(vcf) ]
    samples_ch          // Channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = Channel.empty()

    // now create a tstree file
    TSINFER_MAJOR(
        focal_vcf_ch,
        samples_ch.first()
    )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
