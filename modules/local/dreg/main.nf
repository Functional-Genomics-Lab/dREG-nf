process DREG_RUN {
    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "*dREG*"

    tag "$meta.id"
    memory '50 GB'
    time '48h'
    cpus 4
    accelerator 1

    input:
    tuple val(meta), path(pos_bw), path(neg_bw)

    output:
    tuple val(meta), path("${prefix}.*"), emit: dREG

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    run_dREG.bsh \\
        ${pos_bw} \\
        ${neg_bw} \\
        ${prefix} \\
        ${params.dreg_train} \\
        ${task.cpus} ${$task.accelerator.request}
    """
}
