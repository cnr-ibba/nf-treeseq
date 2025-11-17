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

    // getting focal samples to keep (plink workflow)
    samples_ch = Channel.fromPath( params.sample2fid, checkIfExists: true )

    // call plink subworkflow
    PLINK_EXTRACT(
        ch_samplesheet_plink,
        samples_ch
    )
    ch_versions = ch_versions.mix(PLINK_EXTRACT.out.versions)

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

    if (params.ancestor_method == 'est-sfs') {
        // prepare ancestral samples, call est-sfs and then tsinfer
        // TODO: this option is plink specific for now
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
            PLINK_EXTRACT.out.vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(REFERENCE.out.versions)
    } else if (params.ancestor_method == 'major') {
        // call tsinfer using major alleles as ancestral alleles
        MAJOR(
            PLINK_EXTRACT.out.vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(MAJOR.out.versions)
    } else if (params.ancestor_method == 'custom') {
        // call tsinfer using custom ancestral alleles from user-provided file
        ancestor_ch = Channel.fromPath( params.ancestor_file, checkIfExists: true )

        CUSTOM(
            PLINK_EXTRACT.out.vcf,
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
