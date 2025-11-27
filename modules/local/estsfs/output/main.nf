
process ESTSFS_OUTPUT {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bunop/tskit:fc4ae9345d625e66"

    input:
    tuple val(meta), path(mapping), path(pvalues)

    output:
    tuple val(meta), path("*.ancestral.csv"),     emit: ancestral
    path "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_est_sfs_output \\
        --mapping ${mapping} \\
        --pvalues ${pvalues} \\
        --output ${prefix}.ancestral.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tskitetude: \$(pip show tskitetude | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
