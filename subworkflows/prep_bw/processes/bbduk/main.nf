process BBDUK {
    publishDir "results/bbduk/trim_fastq", mode: 'copy', pattern: "*.trim.fastq.gz"
    publishDir "results/bbduk/trim_stats", mode: 'copy', pattern: "*.trimstats.txt"

    cpus 16
    time '4h'
    memory '32 GB'
    
    input:
    tuple val(meta), path(indices), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: trimmed_reads
    path "*.{refstats,trimstats}.txt", emit: trim_stats

    script:
    prefix_pe = reads[0].toString() - ~/(_1\.)?(_R1)?(\.fq)?(fq)?(\.fastq)?(fastq)?(\.gz)?$/
    prefix_se = reads[0].toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/

    def rnastrandness = ''
    if (params.forwardStranded && !params.unStranded){
        rnastrandness = params.singleEnd ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (params.reverseStranded && !params.unStranded){
        rnastrandness = params.singleEnd ? '--rna-strandness R' : '--rna-strandness RF'
    }

    def bbmap_adapters = params.bbmap_adapters
    def indices_path = params.hisat_indices
    def single_end = meta.single_end
    def name = meta.id

    if (!single_end && params.flip) {
        """
        echo ${name}         
            
        bbduk.sh -Xmx40g \
                t=16 \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}_1.trim.fastq.gz \
                out2=${prefix_pe}_2.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_pe}.trimstats.txt \
                refstats=${prefix_pe}.refstats.txt

        reformat.sh -Xmx40g \
                t=16 \
                in=${prefix_pe}_1.trim.fastq.gz \
                in2=${prefix_pe}_2.trim.fastq.gz \
                out=${prefix_pe}_1.flip.trim.fastq.gz \
                out2=${prefix_pe}_2.flip.trim.fastq.gz \
                rcomp=t                
        """
    } else if (params.singleEnd && params.flip) {
        """
        echo ${name}        
    
        bbduk.sh -Xmx40g \
                t=16 \
                in=${reads} \
                out=${prefix_se}.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
                maq=10 minlen=25 \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_se}.trimstats.txt \
                refstats=${prefix_se}.refstats.txt
                
        reformat.sh -Xmx40g \
                t=16 \
                in=${prefix_se}.trim.fastq.gz \
                out=${prefix_se}.flip.trim.fastq.gz \
                rcomp=t                 
        """
    } else if (!single_end && params.flipR2) {
        """
        echo ${prefix_pe}

        bbduk.sh -Xmx40g \
                t=16 \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}_1.trim.fastq.gz \
                out2=${prefix_pe}_2.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_pe}.trimstats.txt \
                refstats=${prefix_pe}.refstats.txt
            
        reformat.sh -Xmx40g \
                t=16 \
                in=${prefix_pe}_1.trim.fastq.gz \
                in2=${prefix_pe}_2.trim.fastq.gz \
                out=${prefix_pe}_1.flip.trim.fastq.gz \
                out2=${prefix_pe}_2.flip.trim.fastq.gz \
                rcompmate=t          
        """
    } else if (!single_end) {
        """
        echo ${prefix_pe}

        bbduk.sh -Xmx40g \
                t=16 \
                in=${reads[0]} \
                in2=${reads[1]} \
                out=${prefix_pe}_1.trim.fastq.gz \
                out2=${prefix_pe}_2.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_pe}.trimstats.txt \
                refstats=${prefix_pe}.refstats.txt                        
        """
    } else {
        """
        echo ${prefix_se}

        bbduk.sh -Xmx40g \
                t=16 \
                in=${reads} \
                out=${prefix_se}.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=rl trimq=10 k=23 mink=11 hdist=1 \
                maq=10 minlen=25 \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${prefix_se}.trimstats.txt \
                refstats=${prefix_se}.refstats.txt                
        """
    } 
}