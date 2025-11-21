//
// call tsinfer using major alleles as ancestral alleles
//
include { TSINFER_MAJOR } from '../../modules/local/tsinfer/major/main'


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
    ch_versions = ch_versions.mix( TSINFER_MAJOR.out.versions )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
