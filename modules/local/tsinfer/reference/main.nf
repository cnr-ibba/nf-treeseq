
process TSINFER_REFERENCE {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/bunop/tskit:devel"
    containerOptions """${ workflow.containerEngine == 'singularity' ?
        "--bind \${HOME}/.cache/" :
        "--volume \${HOME}/.cache/:/.cache/" }"""

    input:
    tuple val(meta), path(vcf)
    path(sample_file)

    output:
    tuple val(meta), path("*.samples"),     emit: samples
    tuple val(meta), path("*.trees"),       emit: trees
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    create_tstree \\
        --vcf ${vcf} \\
        --focal ${sample_file} \\
        --ancestral_as_reference \\
        --output_samples ${prefix}.samples \\
        --output_trees ${prefix}.trees \\
        --num_threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tskitetude: \$(pip show tskitetude | sed -n 's/^Version: //p')
        tskit: \$(pip show tskit | sed -n 's/^Version: //p')
        tsinfer: \$(pip show tsinfer | sed -n 's/^Version: //p')
        tsdate: \$(pip show tsdate | sed -n 's/^Version: //p')
        tszip: \$(pip show tszip | sed -n 's/^Version: //p')
    END_VERSIONS
    """
}
