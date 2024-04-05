process DREG_PREP {

    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(bam_file), val(index)
    path  sizes
    val assay_type

    output:
    tuple val(meta), path("${prefix}_plus.rpm.bw"), emit: plus_rpm_bigwig
    tuple val(meta), path("${prefix}_minus.rpm.bw"), emit: minus_rpm_bigwig
    tuple val(meta), path("${prefix}_plus.bw"), emit: plus_bigwig
    tuple val(meta), path("${prefix}_minus.bw"), emit: minus_bigwig

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template "proseq2.0"
}
