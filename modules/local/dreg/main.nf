process DREG_RUN {
    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "*dREG*"

    container "docker.io/biohpc/dreg"

    tag "$meta.id"
    memory '50 GB'
    time '48h'
    cpus 16
    accelerator 1

    input:
    tuple val(meta), path(pos_bw), path(neg_bw)
    path model

    output:
    tuple val(meta), path("${prefix}.*"), emit: dREG

    script:
    def prefix = "${meta.id}"
    """
    run_dREG.R \\
        ${pos_bw} \\
        ${neg_bw} \\
        ${prefix} \\
        ${model} \\
        ${task.cpus} \\
        ${task.accelerator.request}
    """
}
