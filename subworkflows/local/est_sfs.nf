//
// call tsinfer by setting ancestrall alleles as defined by est-sfs
//
include { PLINK_SUBSET as ANCIENT_SUBSET    } from '../../modules/local/plink_subset.nf'
include { PLINK_RECODE as ANCIENT_RECODE    } from '../../modules/nf-core/plink/recode/main'
include { BCFTOOLS_NORM as ANCIENT_NORM     } from '../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_SPLIT as ANCIENT_SPLIT   } from '../../modules/nf-core/bcftools/split/main'
include {
    TABIX_TABIX as ANCIENT_TABIX;
    TABIX_TABIX as ANCIENT_SPLIT_TABIX      } from '../../modules/nf-core/tabix/tabix/main'
include { BCFTOOLS_MERGE                    } from '../../modules/nf-core/bcftools/merge/main'
include { ESTSFS_INPUT                      } from '../../modules/local/estsfs_input'
include { ESTSFS                            } from '../../modules/cnr-ibba/estsfs/main'
include { ESTSFS_OUTPUT                     } from '../../modules/local/estsfs_output'
include { TSINFER_ESTSFS                    } from '../../modules/local/tsinfer_estsfs'


process GENERATE_SEED {
    output:
    path 'seedfile.txt'

    '''
    echo $RANDOM > seedfile.txt
    '''
}


workflow EST_SFS {
    take:

    outgroup1           // string: outgroup1 name
    outgroup2           // string: outgroup2 name (can be null)
    outgroup3           // string: outgroup3 name (can be null)
    plink_input_ch      // Channel: plink input files [ meta, Path(bed), Path(bim), Path(fam) ]
    genome_ch           // Channel: genome files [ meta, Path(fasta) ]
    focal_vcf_ch        // Channel: focal vcf file (phased) [ meta, Path(vcf) ]
    focal_tbi_ch        // Channel: focal tbi file (phased) [ meta, Path(tbi) ]
    samples_ch          // Channel: samples file [ meta, Path(samples) ]

    main:

    ch_versions = Channel.empty()

    // collect the outgroup sample list files. At least one outgroup
    outgroup1_ch = Channel.fromPath( outgroup1, checkIfExists: true)
    outgroup2_ch = outgroup2 ? Channel.fromPath(params.outgroup2, checkIfExists: true): Channel.empty()
    outgroup3_ch = outgroup3 ? Channel.fromPath(params.outgroup3, checkIfExists: true): Channel.empty()
    outgroup_files_ch = outgroup1_ch
        .concat(outgroup2_ch)
        .concat(outgroup3_ch)
    outgroups_ch = outgroup_files_ch
        .splitCsv(header: ["breed", "sample_id"], sep: "\t", strip: true)
        .map{ it -> "${it.breed}\t${it.sample_id}" }
        .collectFile(name: 'outgroups.txt', newLine: true)
        // .view()

    ancient_input_ch = plink_input_ch
        .map{ meta, bed, bim, fam -> [[id: "${bed.getBaseName(1)}.ancient"], bed, bim, fam] }
        // .view()

    // extract the ancient samples
    ANCIENT_SUBSET(ancient_input_ch, outgroups_ch)
    ch_versions = ch_versions.mix(ANCIENT_SUBSET.out.versions)

    // transform the plink files to vcf
    ANCIENT_RECODE(ANCIENT_SUBSET.out.bed.join(ANCIENT_SUBSET.out.bim).join(ANCIENT_SUBSET.out.fam))
    ch_versions = ch_versions.mix(ANCIENT_RECODE.out.versions)

    // Normalize ancient VCF
    ANCIENT_NORM(
        // the third element of the input channel is a tbi
        ANCIENT_RECODE.out.vcfgz.map{ meta, vcf -> [meta, vcf, []] },
        genome_ch
    )
    ch_versions = ch_versions.mix(ANCIENT_NORM.out.versions)

    // index ancient vcf
    ANCIENT_TABIX(ANCIENT_NORM.out.vcf)
    ch_versions = ch_versions.mix(ANCIENT_TABIX.out.versions)

    // split data by chromosomes for ancestor
    ANCIENT_SPLIT(ANCIENT_NORM.out.vcf.join(ANCIENT_TABIX.out.tbi))
    ch_versions = ch_versions.mix(ANCIENT_SPLIT.out.versions)

    // index the splitted ancestor
    ANCIENT_SPLIT_TABIX(ANCIENT_SPLIT.out.split_vcf.transpose())
    ch_versions = ch_versions.mix(ANCIENT_SPLIT_TABIX.out.versions)

    // merge the ancient and focal vcf. Prepare input channels
    vcf_ch = focal_vcf_ch
        .map{ meta, it -> [[id: "samples-merged.${meta.chrom}"], it] }
        .concat(
            ANCIENT_SPLIT.out.split_vcf
                .transpose()
                .map{ meta, vcf ->
                    chrom = vcf.name.tokenize(".")[-3]
                    [[id: "samples-merged.${chrom}"], vcf]
                }
        )
        .groupTuple()
        // .view()

    // merge the ancient and focal tbi
    tbi_ch = focal_tbi_ch
        .map{ meta, it -> [[id: "samples-merged.${meta.chrom}"], it] }
        .concat(
            ANCIENT_SPLIT_TABIX.out.tbi
                .map{ meta, tbi ->
                    chrom = tbi.name.tokenize(".")[-4]
                    [[id: "samples-merged.${chrom}"], tbi]
                }
        )
        .groupTuple()
        // .view()

    bcftools_input_ch = vcf_ch.join(tbi_ch)
        // .view()

    // merge the ancient and focal vcf
    BCFTOOLS_MERGE(bcftools_input_ch, [[], []], [[], []], [])
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)

    // calculate ancestral alleles. I need to use the first() method to transform
    // the queue in a value channel
    ESTSFS_INPUT(
        BCFTOOLS_MERGE.out.merged_variants,
        samples_ch.first(),
        outgroup_files_ch.collect()
    )

    // determine a seedfile
    seedfile = GENERATE_SEED()

    // call custom est-sfs
    ESTSFS(
        ESTSFS_INPUT.out.config
            .join(ESTSFS_INPUT.out.input)
            .combine(seedfile)
    )
    ch_versions = ch_versions.mix(ESTSFS.out.versions)

    ESTSFS_OUTPUT(ESTSFS_INPUT.out.mapping.join(ESTSFS.out.pvalues_out))

    tsinfer_in_ch = focal_vcf_ch
        .map{ meta, vcf -> [meta.chrom, meta, vcf] }
        .join(
            ESTSFS_OUTPUT.out.ancestral
                .map{ meta, ancestral ->
                        chrom = ancestral.name.tokenize(".")[-3]
                        [chrom, ancestral]
                },
            by: [0],
            failOnMismatch: true
        ).map{ chrom, meta, vcf, ancestral -> [[id: meta.id], vcf, ancestral]}
        // .view()

    // now create a tstree file
    TSINFER_ESTSFS(
        tsinfer_in_ch,
        samples_ch.first()
    )

    emit:

    versions       = ch_versions                    // channel: [ versions.yml ]
}
