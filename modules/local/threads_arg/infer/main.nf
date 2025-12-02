
process THREADS_INFER {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/bunop/threads_arg:03e31e528ad47bd9"

    input:
    tuple val(meta), path(pgen), path(psam), path(pvar)
    path(demography)

    output:
    tuple val(meta), path("*.threads"), emit: threads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    threads \\
        infer \\
        $args \\
        --num_threads $task.cpus \\
        --pgen $pgen \\
        --demography $demography \\
        --recombination_rate ${params.recombination_rate} \\
        --mutation_rate ${params.mutation_rate} \\
        --query_interval ${params.threads_query_interval} \\
        --fit_to_data \\
        --save_metadata \\
        --out ${prefix}.threads

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

    touch ${prefix}.threads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        threads-arg: \$(pip show threads-arg | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
