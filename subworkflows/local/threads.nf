//
// call threads on VCF data
//
include { BCFTOOLS_VIEW as BCFTOOLS_BIALLELIC   } from '../../modules/nf-core/bcftools/view/main'
include { PLINK2_VCF                            } from '../../modules/nf-core/plink2/vcf/main'
include { GAWK as MAKE_SHAPEIT                  } from '../../modules/nf-core/gawk/main'


workflow THREADS {
    take:
    focal_vcf_ch        // channel: focal vcf file (phased) [ meta, Path(vcf), Path(tbi) ]
    demography_ch       // channel: demography file [ Path(demography) ]

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

    // call gawk to create shapeit format files
    MAKE_SHAPEIT(
        PLINK2_VCF.out.pvar,
        [],
        false
    )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
