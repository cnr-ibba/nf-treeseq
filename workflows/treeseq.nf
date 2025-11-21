/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_nf-treeseq_pipeline'
include { PLINK_EXTRACT                     } from '../subworkflows/local/plink_extract'
include { EST_SFS                           } from '../subworkflows/local/est_sfs'
include { BCFTOOLS_SPLIT as FOCAL_SPLIT     } from '../modules/nf-core/bcftools/split/main'
include { BEAGLE5_BEAGLE as FOCAL_BEAGLE    } from '../modules/nf-core/beagle5/beagle/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_REHEADER                 } from '../modules/nf-core/bcftools/reheader/main'
include { TABIX_TABIX as REHEADER_TABIX     } from '../modules/nf-core/tabix/tabix/main'
include { REFERENCE                         } from '../subworkflows/local/reference'
include { MAJOR                             } from '../subworkflows/local/major'
include { CUSTOM                            } from '../subworkflows/local/custom'

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

    // getting focal samples to keep (this is required to add pupulation information
    // to treesequences output files)
    samples_ch = Channel.fromPath( params.sample2fid, checkIfExists: true )

    // need to define a genome channel
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
        .map{ it -> [[ id: "${it.getBaseName()}" ], it]}
        // .view()

    // this will use VCF received from samplesheet (vcf workflow)
    FOCAL_SPLIT(ch_samplesheet_vcf)
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

    if (params.ancestor_method == 'est-sfs') {
        // prepare ancestral samples, call est-sfs and then tsinfer
        // TODO: this option is plink specific for now
        // call plink subworkflow
        PLINK_EXTRACT(
            ch_samplesheet_plink,
            samples_ch
        )
        ch_versions = ch_versions.mix(PLINK_EXTRACT.out.versions)

        EST_SFS(
            params.outgroup1,
            params.outgroup2,
            params.outgroup3,
            PLINK_EXTRACT.out.plink_input,
            PLINK_EXTRACT.out.genome,
            PLINK_EXTRACT.out.vcf,
            PLINK_EXTRACT.out.tbi,
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
