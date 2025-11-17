/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_nf-treeseq_pipeline'
include { PLINK_SUBSET as FOCAL_SUBSET      } from '../modules/local/plink_subset.nf'
include { PLINK_RECODE as FOCAL_RECODE      } from '../modules/nf-core/plink/recode/main'
include { BCFTOOLS_NORM as FOCAL_NORM       } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SPLIT as FOCAL_SPLIT     } from '../modules/nf-core/bcftools/split/main'
include { BEAGLE5_BEAGLE as FOCAL_BEAGLE    } from '../modules/nf-core/beagle5/beagle/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_REHEADER                 } from '../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX as REHEADER_TABIX     } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE                    } from '../modules/nf-core/bcftools/merge/main'
include { PLINK_EXTRACT                     } from '../subworkflows/local/plink_extract'
include { EST_SFS                           } from '../subworkflows/local/est_sfs'
include { REFERENCE                         } from '../subworkflows/local/reference'
include { MAJOR                             } from '../subworkflows/local/major'
include { CUSTOM                            } from '../subworkflows/local/custom'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TREESEQ {

    take:

    ch_samplesheet_plink // channel: samplesheet for plink files read in from --input
    ch_samplesheet_vcf // channel: samplesheet for vcf files read in from --input

    main:

    ch_versions = Channel.empty()

    // need to define a genome channel
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
        .map{ it -> [[ id: "${it.getBaseName()}" ], it]}
        // .view()

    // getting focal samples to keep (plink workflow)
    samples_ch = Channel.fromPath( params.plink_keep, checkIfExists: true )

    // .view()
    // call plink subworkflow
    PLINK_EXTRACT(
        ch_samplesheet_plink,
        samples_ch,
        genome_ch
    )
    ch_versions = ch_versions.mix(PLINK_EXTRACT.out.versions)

    // split data by chromosomes for focal
    FOCAL_SPLIT(PLINK_EXTRACT.out.vcf.join(PLINK_EXTRACT.out.tbi))
    ch_versions = ch_versions.mix(FOCAL_SPLIT.out.versions)

    // } else if (params.vcf_file) {
    //     // getting input files
    //     vcf_ch = Channel.fromPath( params.vcf_file, checkIfExists: true )
    //         .map{ it -> [[ id: "${it.getBaseName(2)}.focal" ], it] }
    //         // .view()
    //     tbi_ch = Channel.fromPath( params.tbi_file, checkIfExists: true )
    //         .map{ it -> [[ id: "${it.getBaseName(3)}.focal" ], it] }
    //         // .view()

    //     FOCAL_SPLIT(vcf_ch.join(tbi_ch))
    //     ch_versions = ch_versions.mix(FOCAL_SPLIT.out.versions)

    // } else {
    //     error("No valid input file provided")
    // }

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

    if (params.ancestor_method == 'est-sfs') {
        // prepare ancestral samples, call est-sfs and then tsinfer
        // TODO: this option is plink specific for now
        EST_SFS(
            params.outgroup1,
            params.outgroup2,
            params.outgroup3,
            PLINK_EXTRACT.out.plink_input,
            genome_ch,
            BCFTOOLS_REHEADER.out.vcf,
            REHEADER_TABIX.out.tbi,
            samples_ch
        )
        ch_versions = ch_versions.mix(EST_SFS.out.versions)
    } else if (params.ancestor_method == 'reference') {
        // call tsinfer using reference alleles as ancestral alleles
        REFERENCE(
            BCFTOOLS_REHEADER.out.vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(REFERENCE.out.versions)
    } else if (params.ancestor_method == 'major') {
        // call tsinfer using major alleles as ancestral alleles
        MAJOR(
            BCFTOOLS_REHEADER.out.vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(MAJOR.out.versions)
    } else if (params.ancestor_method == 'custom') {
        // call tsinfer using custom ancestral alleles from user-provided file
        ancestor_ch = Channel.fromPath( params.ancestor_file, checkIfExists: true )

        CUSTOM(
            BCFTOOLS_REHEADER.out.vcf,
            samples_ch,
            ancestor_ch
        )
        ch_versions = ch_versions.mix(CUSTOM.out.versions)
    } else {
        error("No valid ancestral allele option provided")
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'nf-treeseq_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    emit:
    versions            = ch_versions                 // channel: [ path(versions.yml) ]
    collated_versions   = ch_collated_versions        // channel: path(collated_versions.yml)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
