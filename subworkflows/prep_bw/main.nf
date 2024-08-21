// include { SRA_DUMP } from './processes/s1_sra_dump/main.nf'
include { BBDUK } from './processes/bbduk/main.nf'
include { MULTI_ALIGN } from './processes/multi_align/main.nf'
include { SAMTOOLS_BAM } from './processes/s3_generate_bam/main.nf'
include { FASTQC } from './processes/fastqc/main.nf'
include { DREG_PREP } from './processes/s5p_dreg_prep/main.nf'


params.hisat_indices = 'hisat2/hg38/genome'

workflow PREP_BW {

    take:
    fastq
    fasta
    fasta_indices
    sizes

    main:
    ch_chrom_sizes = Channel.empty()

    if (!sizes) {
        //
        // Create chromosome sizes file
        //
        CUSTOM_GETCHROMSIZES ( fasta.map { [ [:], it ] } )
        ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    } else {
        ch_chrom_sizes = sizes
    }

    fastq
        | map { [it[0], params.hisat_indices, it[1]] }
        | BBDUK // -> { trimmed_reads, trim_stats }
        | set { ch_fastq_trimmed }

    ch_fastq_trimmed.trimmed_reads
        | map { it[1] } // strip off name and just pass the reads
        | FASTQC // -> { fastqc_results }
    
    Channel.fromList([tuple(aligner: 'bwa-mem', index: fasta), tuple(aligner: 'bowtie2', index: fasta), tuple(aligner: 'hisat2', index: params.hisat_indices)])
        | set { aligners }

    ch_fastq_trimmed.trimmed_reads
        | combine(aligners)
        | map { [
                    it[0] + [ aligner : it[2].aligner ],
                    it[2].aligner, 
                    it[1], 
                    (it[2].index ? it[2].index : fasta), 
                    (it[2].index ? "index" : "genome")
                ] }
        | set { multi_align_input }
    
    MULTI_ALIGN(multi_align_input, Channel.fromPath(fasta_indices).collect())
        | set { ch_multi_align_sam }


    ch_multi_align_sam.aligned_sam
        | map { [it[0], fasta, it[1]] } // -> { meta, genome_path, mapped_sam }
        | SAMTOOLS_BAM  // -> { sorted_bam, sorted_bam_index, sorted_bam_flagstat, sorted_cram, sorted_cram_index }
        | set { ch_bam_and_cram }

    ch_bam_and_cram.sorted_cram.
        merge(SAMTOOLS_BAM.out.sorted_cram_index)
        | map { [it[0], it[1], it[3]] } // -> tuple val(meta), file(cram), file(cram_idx)
        | set { dreg_prep_input }

    DREG_PREP(dreg_prep_input, sizes)

    emit:
    DREG_PREP.out.plus_bigwig
    | join (DREG_PREP.out.minus_bigwig, by: [0])
}
