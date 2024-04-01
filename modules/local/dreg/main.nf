process dreg_run {
    println "Log[6]: Running dREG"
    println "Log[6]: N.B. Requires GPUs"

    tag "$prefix"
    memory '50 GB'
    time '48h'
    cpus 4
    queue 'titan'
    clusterOptions '--gres=gpu'

    publishDir "${params.outdir}/dreg/", mode: 'copy', pattern: "*dREG*"

    when:
    params.dreg

    input:
    tuple val(prefix), file(pos_bw), file(neg_bw) from dreg_bigwig

    output:
    tuple val(prefix), file ("${prefix}.*") into dREG_out

    script:
        """
        bash ${params.dreg_path} \
	     ${pos_bw} \
	     ${neg_bw} \
	     ${prefix} \
	     ${params.dreg_train} \
	     4 1
        """
}
