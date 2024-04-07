process PROSEQ2 {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)
    path bwa_index
    path chromInfo
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
    def reads_command = meta.single_end ? "-SE" : "-PE"
    def required_se_options = meta.single_end ? assay_type == "GROseq" ? "-G" : "-P" : ""
    // TODO PE
    """
    mv $reads ${prefix}.fastq.gz

    proseq2.0.bsh \\
        $reads_command \\
        -i $bwa_index \\
        -c $chromInfo \\
        -I $prefix \\
        $required_se_options \\
        -4DREG \\
        --thread=${task.cpus}
    """
}
