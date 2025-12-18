//
// call tsinfer using major alleles as ancestral alleles
//
include { TSINFER_MAJOR } from '../../modules/local/tsinfer/major/main'


workflow MAJOR {
    take:
    focal_vcf_ch        // channel: focal vcf file (phased) [ meta, Path(vcf) ]
    samples_ch          // channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = channel.empty()

    // now create a tstree file
    TSINFER_MAJOR(
        focal_vcf_ch,
        samples_ch.first()
    )
    ch_versions = ch_versions.mix( TSINFER_MAJOR.out.versions )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
