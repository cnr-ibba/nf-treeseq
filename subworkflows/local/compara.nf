//
// call tsinfer using reference alleles as ancestral alleles
//
include { TSINFER_COMPARA } from '../../modules/local/tsinfer_compara'


workflow COMPARA {
    take:
    focal_vcf_ch        // Channel: focal vcf file (phased) [ meta, Path(vcf) ]
    samples_ch          // Channel: samples file [ Path(samples) ]
    ancestor_ch         // Channel: ancestral file [ Path(ancestor) ]

    main:

    ch_versions = Channel.empty()

    // now create a tstree file
    TSINFER_COMPARA(
        focal_vcf_ch,
        samples_ch.first(),
        ancestor_ch.first()
    )

    emit:

    versions       = ch_versions                    // channel: [ versions.yml ]
}
