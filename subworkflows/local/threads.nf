//
// call threads on VCF data
//
include { BCFTOOLS_VIEW as BCFTOOLS_BIALLELIC   } from '../../modules/nf-core/bcftools/view/main'
include { PLINK2_VCF                            } from '../../modules/nf-core/plink2/vcf/main'
include { THREADS_INFER                         } from '../../modules/local/threads_arg/infer/main'
include { THREADS_CONVERT                       } from '../../modules/local/threads_arg/convert/main'


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
    ch_versions = ch_versions.mix( PLINK2_VCF.out.versions )

    // running threads
    THREADS_INFER(
        PLINK2_VCF.out.pgen
            .join(PLINK2_VCF.out.psam)
            .join(PLINK2_VCF.out.pvar),
        demography_ch.first()
    )
    ch_versions = ch_versions.mix( THREADS_INFER.out.versions )

    // packing outputs in tree sequence format (compressed)
    THREADS_CONVERT(
        THREADS_INFER.out.threads
    )
    ch_versions = ch_versions.mix( THREADS_CONVERT.out.versions )

    emit:
    versions       = ch_versions                    // channel: [ versions.yml ]
}
