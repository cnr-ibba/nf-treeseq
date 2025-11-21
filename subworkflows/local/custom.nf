//
// call tsinfer using reference alleles as ancestral alleles
//
include { TSINFER_CUSTOM } from '../../modules/local/tsinfer/custom/main'


workflow CUSTOM {
    take:
    focal_vcf_ch        // channel: focal vcf file (phased) [ meta, Path(vcf) ]
    samples_ch          // channel: samples file [ Path(samples) ]
    ancestor_ch         // channel: ancestral file [ Path(ancestor) ]

    main:

    ch_versions = channel.empty()

    // now create a tstree file
    TSINFER_CUSTOM(
        focal_vcf_ch,
        samples_ch.first(),
        ancestor_ch.first()
    )
    ch_versions = ch_versions.mix( TSINFER_CUSTOM.out.versions )

    emit:

    versions       = ch_versions                    // channel: [ versions.yml ]
}
