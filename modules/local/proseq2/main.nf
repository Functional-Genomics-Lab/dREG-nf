process PROSEQ2 {
    tag "$meta.id"
    label 'process_high'
    label 'process_long'

    errorStrategy 'finish'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)
    path bwa_index
    path chromInfo
    val assay_type

    output:
    tuple val(meta), path("results/${prefix}_QC_plus.rpm.bw"), emit: plus_rpm_bigwig, optional: true
    tuple val(meta), path("results/${prefix}_QC_minus.rpm.bw"), emit: minus_rpm_bigwig, optional: true
    tuple val(meta), path("results/${prefix}_QC_plus.bw"), emit: plus_bigwig
    tuple val(meta), path("results/${prefix}_QC_minus.bw"), emit: minus_bigwig

    when:
    task.ext.when == null || task.ext.when

    script:
    mata.aligner = "proseq2.0"
    prefix = task.ext.prefix ?: "${meta.id}.${meta.aligner}"
    def reads_command = meta.single_end ? "-SE" : "-PE"
    def required_se_options = meta.single_end ? assay_type == "GROseq" ? "-G" : "-P" : ""

    // TODO PE
    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

    mv $reads ${prefix}.fastq.gz
    mkdir tmp

    proseq2.0.bsh \\
        $reads_command \\
        -i \$INDEX \\
        -c $chromInfo \\
        -I $prefix \\
        $required_se_options \\
        -4DREG \\
        --thread=${task.cpus} \\
        -T ./tmp \\
        -O results
    """
}
