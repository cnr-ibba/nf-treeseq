//
// Extract and phase focal samples from a VCF file
//
include { BCFTOOLS_SPLIT as FOCAL_SPLIT     } from '../../modules/nf-core/bcftools/split/main'
include { BEAGLE5_BEAGLE as FOCAL_BEAGLE    } from '../../modules/nf-core/beagle5/beagle/main'
include { SAMTOOLS_FAIDX                    } from '../../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_REHEADER                 } from '../../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX as REHEADER_TABIX     } from '../../modules/nf-core/tabix/tabix/main'

workflow VCF_EXTRACT {
    take:
    ch_samplesheet          // Channel: vcf input prefix [ meta, vcf, index ]
    samples_ch              // Channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = Channel.empty()

    // need to define a genome channel
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
        .map{ it -> [[ id: "${it.getBaseName()}" ], it]}
        // .view()

    // split VCF by chromosomes
    FOCAL_SPLIT(ch_samplesheet)
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
    genome          = genome_ch
    versions        = ch_versions                    // channel: [ versions.yml ]
}
