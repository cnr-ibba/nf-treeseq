//
// call threads on VCF data
//
include { BCFTOOLS_VIEW as BCFTOOLS_BIALLELIC   } from '../../modules/nf-core/bcftools/view/main'
include { PLINK2_VCF                            } from '../../modules/nf-core/plink2/vcf/main'


workflow THREADS {
    take:
    focal_vcf_ch        // channel: focal vcf file (phased) [ meta, Path(vcf), Path(tbi) ]

    main:

    ch_versions = channel.empty()

    // keep only biallelic sites
    BCFTOOLS_BIALLELIC(
        focal_vcf_ch,
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix( BCFTOOLS_BIALLELIC.out.versions )

    // call plink2 to create pgen/psam/pvar files
    PLINK2_VCF(
        BCFTOOLS_BIALLELIC.out.vcf
    )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
