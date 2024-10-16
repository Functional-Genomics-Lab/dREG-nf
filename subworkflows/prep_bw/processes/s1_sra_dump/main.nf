process SRA_DUMP {
    publishDir "results/s1_sra_dump", mode: 'rellink'

    input:
    file input_sra

    output:
    path '*.fastq.gz', emit: fastq

    script:
    prefix = input_sra.baseName
    if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}
        fastq-dump --split-3 ${reads} --gzip
        """
    } else if (!params.threadfqdump) {
        """
        echo ${prefix}
        fastq-dump ${reads} --gzip
        """
    } else if (!params.singleEnd) {
        """
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --split-3 \
            --sra-id ${reads}
        """
    } else {
        """
        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --sra-id ${reads}
        """
    }
}