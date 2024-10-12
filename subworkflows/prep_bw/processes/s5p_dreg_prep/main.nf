process DREG_PREP {
    publishDir "results/s5p_dreg_prep", mode: 'rellink'

    cpus 1
    memory '50 GB'

    input:
    tuple val(meta), path(cram_file), path(cram_idx)
    path(chrom_sizes)

    output:
    tuple val(meta), path("*.pos.bw"), emit: plus_bigwig
    tuple val(meta), path("*.neg.bw"), emit: minus_bigwig

    script:
    def genome = params.fasta
    def name = "${meta.id}.${meta.aligner}"
    """
    echo "Creating BigWigs suitable as inputs to dREG"
    
    export CRAM_REFERENCE=${genome}    
    
    bamToBed -i ${cram_file} | awk 'BEGIN{OFS="\t"} (\$5 > 0){print \$0}' | \
    awk 'BEGIN{OFS="\t"} (\$6 == "+") {print \$1,\$2,\$2+1,\$4,\$5,\$6}; (\$6 == "-") {print \$1, \$3-1,\$3,\$4,\$5,\$6}' \
    > ${name}.dreg.bed
    sortBed -i ${name}.dreg.bed > ${name}.dreg.sort.bed
    
    echo positive strand processed to bedGraph
    
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand + \
            > ${name}.pos.bedGraph
    sortBed \
            -i ${name}.pos.bedGraph \
            > ${name}.pos.sort.bedGraph
            
    bedtools genomecov \
            -bg \
            -i ${name}.dreg.sort.bed \
            -g ${chrom_sizes} \
            -strand - \
            | awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' > ${name}.neg.bedGraph
    sortBed \
            -i ${name}.neg.bedGraph \
            > ${name}.neg.sort.bedGraph
            
    echo negative strand processed to bedGraph

    for strand in pos neg; do 
        if [ -s "${name}.\${strand}.sort.bedGraph" ]; then
            echo "Creating \${strand} strand bigwig"
            bedGraphToBigWig "${name}.\${strand}.sort.bedGraph" "${chrom_sizes}" "${name}.\${strand}.bw"
        else
            echo "\${strand} strand bedGraph is empty. Creating stub bigwig"
            touch "${name}.\${strand}.bw"
            truncate -s 0 "${name}.\${strand}.bw"
        fi
    done
    
    echo bedGraph to bigwig done
    """
}