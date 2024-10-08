
process TSINFER_ESTSFS {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/bunop/tskit:0.1.2"
    containerOptions """${ workflow.containerEngine == 'singularity' ?
        "--bind ${HOME}/.cache/" :
        "--volume ${HOME}/.cache/:/.cache/" }"""

    input:
    tuple val(meta), path(vcf), path(ancestral)
    path(sample_file)

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
        --ancestral_estsfs ${ancestral} \\
        --output_samples ${prefix}.samples \\
        --output_trees ${prefix}.trees \\
        --num_threads $task.cpus \\
        $args
    """
}
