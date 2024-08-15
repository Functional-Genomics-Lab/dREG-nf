process DREG_RUN {
    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "dREG_output_*"

    container "${moduleDir}/choose_dreg_image.sh".execute().text.trim()

    tag "$meta.id"
    memory '50 GB'
    time '48h'
    cpus 16
    accelerator 1

    input:
    tuple val(meta), path(pos_bw), path(neg_bw)
    path model

    output:
    tuple val(meta), path("dREG_output_*"), emit: dREG

    script:
    def prefix = "dREG_output_${meta.id}_"
    """
    /dREG/run_dREG.R \\
        ${pos_bw} \\
        ${neg_bw} \\
        ${prefix} \\
        ${model} \\
        ${task.cpus} \\
        0
    """
}
