process FASTQC {
    publishDir "results/fastqc", mode: 'symlink'
    
    cpus 1
    memory '8 GB'

    input:
    file input_fastq

    output:
    tuple val(input_fastq.baseName), file("*.html"), file("*.zip"), emit: fastqc_results

    script:
    prefix = input_fastq.baseName
    """
    echo ${prefix}
    fastqc ${input_fastq}
    """
}