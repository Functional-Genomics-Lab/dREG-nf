process SAMTOOLS_BAM {
    publishDir "results/s3_generate_bam", mode: 'rellink'

    cpus 16
    memory '32 GB'

    input:
    tuple val(meta), val(genome_path), file(mapped_sam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sorted_bam
    tuple val(meta), path("*.sorted.bam.bai"), emit: sorted_bam_index
    tuple val(meta), path("*.flagstat"), emit: sorted_bam_flagstat
    tuple val(meta), path("*.sorted.cram"), emit: sorted_cram
    tuple val(meta), path("*.sorted.cram.crai"), emit: sorted_cram_index

    script:
    def name = meta.id
    """
    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x904 -c ${name}.sorted.bam > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome_path} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
}