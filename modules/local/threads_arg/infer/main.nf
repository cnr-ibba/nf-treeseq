
process THREADS_INFER {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/bunop/threads_arg:894177275e63d259"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar), path(map)
    path(demography)

    output:
    tuple val(meta), path("*.threads"), emit: threads
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
        infer \\
        $args \\
        --num_threads $task.cpus \\
        --pgen $pgen \\
        --map $map \\
        --demography $demography \\
        --mutation_rate ${params.threads_mutation_rate} \\
        --query_interval ${params.threads_query_interval} \\
        --out ${prefix}.threads

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

    touch ${prefix}.threads
    """
}
