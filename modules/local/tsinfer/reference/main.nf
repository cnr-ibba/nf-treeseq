
process TSINFER_REFERENCE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_long'

    container "docker.io/bunop/tskit:0.5.1"
    containerOptions """${ workflow.containerEngine == 'singularity' ?
        "--bind \${HOME}/.cache/" :
        "--volume \${HOME}/.cache/:/.cache/" }"""

    input:
    tuple val(meta), path(vcf)
    path(sample_file)

    output:
    tuple val(meta), path("*.samples"),     emit: samples
    tuple val(meta), path("*.trees.tsz"),   emit: trees
    path "versions.yml",                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    mkfifo ${prefix}.trees

    # this is a cleanup function to ensure that background processes are killed
    # and temporary files are removed
    cleanup() {
        kill \$CREATE_PID \$TSZIP_PID 2>/dev/null || true
        rm -f ${prefix}.trees
    }
    trap cleanup EXIT INT TERM

    # call create_tstree and direct its output to a named pipe
    create_tstree \\
        --vcf ${vcf} \\
        --focal ${sample_file} \\
        --ancestral_as_reference \\
        --output_samples ${prefix}.samples \\
        --output_trees ${prefix}.trees \\
        --num_threads $task.cpus \\
        $args &

    # collect the PID of the background process
    CREATE_PID=\$!

    # compress the trees output from the named pipe and collect its PID
    tszip ${prefix}.trees &
    TSZIP_PID=\$!

    # Wait for create_tstree
    wait \$CREATE_PID
    CREATE_EXIT=\$?

    # Check if create_tstree succeeded
    if [ \$CREATE_EXIT -ne 0 ]; then
        echo "ERROR: create_tstree failed" >&2
        exit \$CREATE_EXIT
    fi

    # wait for tszip to finish
    wait \$TSZIP_PID

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
