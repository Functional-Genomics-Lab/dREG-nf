process MULTI_ALIGN {
    cpus 16
    time '4h'
    memory '32 GB'

    input:
    tuple val(meta), val(aligner), path(reads), val(reference), val(reference_type)
    path indices

    output:
    tuple val(meta), path("*.sam"), emit: aligned_sam

    script:
    def forward_opt = "--input-forward ${reads[0]}"
    def reverse_opt = reads.size() == 2 ? "--input-reverse ${reads[1]}" : ""
    def aligner = aligner.toLowerCase()
    def input_flag = reference_type == "genome" ? "--genome ${reference}" : "--index ${reference}"
    
    """
    multi_align \
        ${forward_opt} \
        ${reverse_opt} \
        --output ${meta.id}.${aligner}.aligned.sam \
        --threads 16 \
        --allow-spliced-alignments false \
        ${aligner} \
        ${input_flag}
    """
}