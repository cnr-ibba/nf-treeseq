
process THREADS_CONVERT {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bunop/threads_arg:03e31e528ad47bd9"

    input:
    tuple val(meta), path(threads)

    output:
    tuple val(meta), path("*.tsz"), emit: tree
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    threads \\
        convert \\
        $args \\
        --add_mutations \\
        --threads ${threads} \\
        --tsz ${prefix}.tsz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        threads-arg: \$(pip show threads-arg | sed -n 's/^Version: //p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tsz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        threads-arg: \$(pip show threads-arg | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
