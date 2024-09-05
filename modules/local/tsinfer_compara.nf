
process TSINFER_COMPARA {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/bunop/tskit:0.1.2"
    containerOptions """${ workflow.containerEngine == 'singularity' ?
        "--bind ${HOME}/.cache/" :
        "--volume ${HOME}/.cache/:/.cache/" }"""

    input:
    tuple val(meta), path(vcf)
    path(sample_file)
    path(ancestor_file)

    output:
    tuple val(meta), path("*.samples"),     emit: samples
    tuple val(meta), path("*.trees"),       emit: trees

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    create_tstree \\
        --vcf ${vcf} \\
        --focal ${sample_file} \\
        --ancestral_ensembl ${ancestor_file} \\
        --output_samples ${prefix}.samples \\
        --output_trees ${prefix}.trees \\
        --num_threads $task.cpus \\
        $args
    """
}
