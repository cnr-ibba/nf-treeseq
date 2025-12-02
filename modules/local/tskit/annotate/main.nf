
process TSKIT_ANNOTATE {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bunop/tskit:7fe0ea482d8c3d98"
    containerOptions """${ workflow.containerEngine == 'singularity' ?
        "--bind \${HOME}/.cache/" :
        "--volume \${HOME}/.cache/:/.cache/" }"""

    input:
    tuple val(meta), path(tree)
    path(sample_file)
    tuple val(software_name), val(software_version)

    output:
    tuple val(meta), path("*.trees.tsz"),   emit: trees
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.trees"
    """
    annotate_tree.py \\
        --input_tsz ${tree} \\
        --sample_file ${sample_file} \\
        --output_tsz ${prefix}.tsz \\
        --software_name ${software_name} \\
        --software_version ${software_version} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tskit: \$(pip show tskit | sed -n 's/^Version: //p')
        tszip: \$(pip show tszip | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.trees"
    """
    echo $args

    touch ${prefix}.tsz
    touch versions.yml
    """
}
