//
// extract samples from plink files and create a indexed VCF
//
include { PLINK_SUBSET                      } from '../../modules/local/plink_subset.nf'
include { PLINK_RECODE                      } from '../../modules/nf-core/plink/recode/main'
include { BCFTOOLS_NORM                     } from '../../modules/nf-core/bcftools/norm/main'
include { TABIX_TABIX as BCFTOOLS_TABIX     } from '../../modules/nf-core/tabix/tabix/main'
include { TABIX_TABIX as REHEADER_TABIX     } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_SPLIT as FOCAL_SPLIT     } from '../../modules/nf-core/bcftools/split/main'
include { BEAGLE5_BEAGLE as FOCAL_BEAGLE    } from '../../modules/nf-core/beagle5/beagle/main'
include { SAMTOOLS_FAIDX                    } from '../../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_REHEADER                 } from '../../modules/nf-core/bcftools/reheader/main'

workflow PLINK_EXTRACT {
    take:

    ch_samplesheet          // Channel: plink input prefix [ meta, prefix ]
    samples_ch              // Channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = Channel.empty()

    // need to define a genome channel
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
        .map{ it -> [[ id: "${it.getBaseName()}" ], it]}
        // .view()

    // here are samples in input
    ch_samplesheet
        .map{ meta, plink_prefix ->
            def bed = file("${plink_prefix}.bed")
            def bim = file("${plink_prefix}.bim")
            def fam = file("${plink_prefix}.fam")
            if ( !bed.exists() || !bim.exists() || !fam.exists() ) {
                error("PLINK binary files (.bed, .bim, .fam) for sample '${meta.id}' not found at: ${plink_prefix}.bed, ${plink_prefix}.bim, ${plink_prefix}.fam")
            }
            return [ meta, [ bed, bim, fam ] ]
        }
        .map{ _meta, plink -> [[ id: "${plink[0].getBaseName(1)}.focal" ], plink[0], plink[1], plink[2]] }
        .set { plink_input_ch }

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
    BCFTOOLS_TABIX(BCFTOOLS_NORM.out.vcf)
    ch_versions = ch_versions.mix(BCFTOOLS_TABIX.out.versions)

    // split data by chromosomes for focal
    FOCAL_SPLIT(BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_TABIX.out.tbi))
    ch_versions = ch_versions.mix(FOCAL_SPLIT.out.versions)

    // get the chromosome name from the vcf file name
    beagle_in_ch = FOCAL_SPLIT.out.split_vcf
        .transpose()
        .map{ meta, vcf ->
            def chrom = vcf.name.tokenize(".")[-3]
            [[id: "${meta.id}.${chrom}", chrom: chrom], vcf, [], [], [], [], [], []]
        }
        // .view()

    // phase and impute with beagle5
    FOCAL_BEAGLE(beagle_in_ch)
    ch_versions = ch_versions.mix(FOCAL_BEAGLE.out.versions)

    // index genome sequence
    SAMTOOLS_FAIDX(genome_ch, [[], []], [])
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // when I have two queues of different size, I can use the first() method
    // to transform the queue in a value channel
    // https://training.nextflow.io/basic_training/channels/#value-channels
    BCFTOOLS_REHEADER(
        FOCAL_BEAGLE.out.vcf.map{ meta, vcf -> [meta, vcf, [], []] },
        SAMTOOLS_FAIDX.out.fai.first()
    )
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions)

    // index beagle genotype
    REHEADER_TABIX(BCFTOOLS_REHEADER.out.vcf)
    ch_versions = ch_versions.mix(REHEADER_TABIX.out.versions)

    emit:

    vcf             = BCFTOOLS_REHEADER.out.vcf
    tbi             = REHEADER_TABIX.out.tbi
    plink_input     = plink_input_ch
    genome          = genome_ch
    versions        = ch_versions                    // channel: [ versions.yml ]
}
