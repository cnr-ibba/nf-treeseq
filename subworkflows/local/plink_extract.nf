//
// extract samples from plink files and create a indexed VCF
//
include { PLINK_SUBSET  } from '../../modules/local/plink_subset.nf'
include { PLINK_RECODE  } from '../../modules/nf-core/plink/recode/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX   } from '../../modules/nf-core/tabix/tabix/main'

workflow PLINK_EXTRACT {
    take:

    plink_input_ch      // Channel: plink input files [ meta, Path(bed), Path(bim), Path(fam) ]
    samples_ch          // Channel: samples file [ meta, Path(samples) ]
    genome_ch           // Channel: genome files [ meta, Path(fasta) ]

    main:

    ch_versions = Channel.empty()

    // extract the samples I want. See modules.config for other options
    PLINK_SUBSET(plink_input_ch, samples_ch)
    ch_versions = ch_versions.mix(PLINK_SUBSET.out.versions)

    // transform the plink files to vcf
    PLINK_RECODE(PLINK_SUBSET.out.bed.join(PLINK_SUBSET.out.bim).join(PLINK_SUBSET.out.fam))
    ch_versions = ch_versions.mix(PLINK_RECODE.out.versions)

    // Normalize focal VCF
    BCFTOOLS_NORM(
        PLINK_RECODE.out.vcfgz.map{ meta, vcf -> [meta, vcf, []] },
        genome_ch
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions)

    // index focal vcf
    TABIX_TABIX(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    emit:

    vcf             = BCFTOOLS_NORM.out.vcf
    tbi             = TABIX_TABIX.out.tbi
    versions        = ch_versions                    // channel: [ versions.yml ]
}
