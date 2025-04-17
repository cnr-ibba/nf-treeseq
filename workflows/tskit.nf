/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowTskit.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { PLINK_SUBSET as FOCAL_SUBSET      } from '../modules/local/plink_subset.nf'
include { PLINK_RECODE as FOCAL_RECODE      } from '../modules/nf-core/plink/recode/main'
include { BCFTOOLS_NORM as FOCAL_NORM       } from '../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SPLIT as FOCAL_SPLIT     } from '../modules/nf-core/bcftools/split/main'
include { BEAGLE5_BEAGLE as FOCAL_BEAGLE    } from '../modules/nf-core/beagle5/beagle/main'
include { SAMTOOLS_FAIDX                    } from '../modules/nf-core/samtools/faidx/main'
include { BCFTOOLS_REHEADER                 } from '../modules/nf-core/bcftools/reheader/main'
include {
    TABIX_TABIX as FOCAL_TABIX;
    TABIX_TABIX as REHEADER_TABIX           } from '../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE                    } from '../modules/nf-core/bcftools/merge/main'
include { EST_SFS                           } from '../subworkflows/local/est_sfs'
include { REFERENCE                         } from '../subworkflows/local/reference'
include { MAJOR                             } from '../subworkflows/local/major'
include { COMPARA                           } from '../subworkflows/local/compara'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow TSKIT {
    ch_versions = Channel.empty()

    // getting plink input files
    bed =  Channel.fromPath( "${params.plink_bfile}.bed" )
    bim =  Channel.fromPath( "${params.plink_bfile}.bim" )
    fam =  Channel.fromPath( "${params.plink_bfile}.fam" )

    plink_input_ch = bed.concat(bim, fam)
        .collect()
        .map{ it -> [[ id: "${it[0].getBaseName(1)}.focal" ], it[0], it[1], it[2]] }
        // .view()

    // getting focal samples to keep
    samples_ch = Channel.fromPath( params.plink_keep, checkIfExists: true )

    // extract the samples I want. See modules.confing for other options
    FOCAL_SUBSET(plink_input_ch, samples_ch)
    ch_versions = ch_versions.mix(FOCAL_SUBSET.out.versions)

    // transform the plink files to vcf
    FOCAL_RECODE(FOCAL_SUBSET.out.bed.join(FOCAL_SUBSET.out.bim).join(FOCAL_SUBSET.out.fam))
    ch_versions = ch_versions.mix(FOCAL_RECODE.out.versions)

    // need to define a genome channel
    genome_ch = Channel.fromPath(params.genome, checkIfExists: true)
        .map{ it -> [[ id: "${it.getBaseName()}" ], it]}
        // .view()

    // Normalize focal VCF
    FOCAL_NORM(
        FOCAL_RECODE.out.vcfgz.map{ meta, vcf -> [meta, vcf, []] },
        genome_ch
    )
    ch_versions = ch_versions.mix(FOCAL_NORM.out.versions)

    // index focal vcf
    FOCAL_TABIX(FOCAL_NORM.out.vcf)
    ch_versions = ch_versions.mix(FOCAL_TABIX.out.versions)

    // split data by chromosomes for focal
    FOCAL_SPLIT(FOCAL_NORM.out.vcf.join(FOCAL_TABIX.out.tbi))
    ch_versions = ch_versions.mix(FOCAL_SPLIT.out.versions)

    // get the chromosome name from the vcf file name
    beagle_in_ch = FOCAL_SPLIT.out.split_vcf
        .transpose()
        .map{ meta, vcf ->
            def chrom = vcf.name.tokenize(".")[-3]
            [[id: "${meta.id}.${chrom}", chrom: chrom], vcf]
        }
        // .view()

    // phase and inpute with beagle5
    FOCAL_BEAGLE(beagle_in_ch, [], [], [], [])
    ch_versions = ch_versions.mix(FOCAL_BEAGLE.out.versions)

    // index genome sequence
    SAMTOOLS_FAIDX(genome_ch, [[], []])
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
    ch_versions = ch_versions.mix(FOCAL_TABIX.out.versions)

    if (params.with_estsfs) {
        // prepare ancestral samples, call est-sfs and then tsinfer
        EST_SFS(
            params.outgroup1,
            params.outgroup2,
            params.outgroup3,
            plink_input_ch,
            genome_ch,
            BCFTOOLS_REHEADER.out.vcf,
            REHEADER_TABIX.out.tbi,
            samples_ch
        )
        ch_versions = ch_versions.mix(EST_SFS.out.versions)
    } else if (params.reference_ancestor) {
        // call tsinfer using reference alleles as ancestral alleles
        REFERENCE(
            BCFTOOLS_REHEADER.out.vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(REFERENCE.out.versions)
    } else if (params.reference_major) {
        // call tsinfer using major alleles as ancestral alleles
        MAJOR(
            BCFTOOLS_REHEADER.out.vcf,
            samples_ch
        )
        ch_versions = ch_versions.mix(MAJOR.out.versions)
    } else if (params.compara_ancestor) {
        // call tsinfer using ancestral alleles from ensembl-compara
        ancestor_ch = Channel.fromPath( params.compara_ancestor, checkIfExists: true )

        COMPARA(
            BCFTOOLS_REHEADER.out.vcf,
            samples_ch,
            ancestor_ch
        )
        ch_versions = ch_versions.mix(COMPARA.out.versions)
    } else {
        error("No valid ancestral allele option provided")
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
