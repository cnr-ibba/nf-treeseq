//
// call tsinfer using reference alleles as ancestral alleles
//

workflow REFERENCE {
    take:

    main:

    ch_versions = Channel.empty()

    emit:

    versions       = ch_versions                    // channel: [ versions.yml ]
}
