/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                  } from 'plugin/nf-schema'
include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText            } from '../subworkflows/local/utils_nfcore_nf-treeseq_pipeline'
include { PLINK_EXTRACT                     } from '../subworkflows/local/plink_extract'
include { VCF_EXTRACT                       } from '../subworkflows/local/vcf_extract'
include { EST_SFS                           } from '../subworkflows/local/est_sfs'
include { REFERENCE                         } from '../subworkflows/local/reference'
include { MAJOR                             } from '../subworkflows/local/major'
include { CUSTOM                            } from '../subworkflows/local/custom'
include { THREADS                           } from '../subworkflows/local/threads'

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

    ch_versions = channel.empty()

    // getting focal samples to keep (this is required to add pupulation information
    // to treesequences output files)
    samples_ch = channel.fromPath( params.sample2fid, checkIfExists: true )

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
    } else {
        // prepare ancestral samples and then tsinfer
        VCF_EXTRACT(
            ch_samplesheet_vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(VCF_EXTRACT.out.versions)

        if (params.ancestor_method == 'reference') {
            // call tsinfer using reference alleles as ancestral alleles
            REFERENCE(
                VCF_EXTRACT.out.vcf,
                samples_ch
            )
            ch_versions = ch_versions.mix(REFERENCE.out.versions)
        } else if (params.ancestor_method == 'major') {
            // call tsinfer using major alleles as ancestral alleles
            MAJOR(
                VCF_EXTRACT.out.vcf,
                samples_ch
            )
            ch_versions = ch_versions.mix(MAJOR.out.versions)
        } else if (params.ancestor_method == 'custom') {
            // call tsinfer using custom ancestral alleles from user-provided file
            ancestor_ch = channel.fromPath( params.ancestor_file, checkIfExists: true )

            CUSTOM(
                VCF_EXTRACT.out.vcf,
                samples_ch,
                ancestor_ch
            )
            ch_versions = ch_versions.mix(CUSTOM.out.versions)
        } else if (params.ancestor_method == 'threads') {
            // create tree sequences using https://palamaralab.github.io/software/threads/
            // open demography file (if provided) or create a default one using Ne
            if (params.threads_demography_file) {
                demography_ch = channel.fromPath(
                    params.threads_demography_file,
                    checkIfExists: true
                )
            } else {
                // create a default demography file using Ne
                demography_ch = channel.of("0 ${params.threads_ne}")
                    | map { content ->
                        def tsv_file = file("${workDir}/demography_default.tsv")
                        tsv_file.text = content
                        return tsv_file
                    }
            }

            // call threads subworkflow
            THREADS(
                VCF_EXTRACT.out.vcf.join(VCF_EXTRACT.out.tbi),
                demography_ch,
                samples_ch
            )
            ch_versions = ch_versions.mix(THREADS.out.versions)
        } else {
            error("No valid ancestral allele option provided")
        }
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
