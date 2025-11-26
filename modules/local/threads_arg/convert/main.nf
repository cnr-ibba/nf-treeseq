
process THREADS_CONVERT {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/bunop/threads_arg:03e31e528ad47bd9"

    input:
    tuple val(meta), path(threads)

    output:
    tuple val(meta), path("*.tsz"), emit: threads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = "0.2.1"
    """
    threads \\
        convert \\
        $args \\
        --threads ${threads} \\
        --tsz ${prefix}.tsz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        threads: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tsz
    """
}
