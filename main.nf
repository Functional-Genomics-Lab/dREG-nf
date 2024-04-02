#!/usr/bin/env nextflow

params.bams = null
params.sizes = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.dict"
params.fasta = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.assay_type = "GROseq"
params.dreg_model = "https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata"

include { CUSTOM_GETCHROMSIZES } from './modules/nf-core/custom/getchromsizes/main'
include { SAMTOOLS_SORT } from './modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'

include { DREG_PREP } from './modules/local/dreg_prep/main'
include { DREG_RUN } from './modules/local/dreg/main'

workflow {
    ch_chrom_sizes = Channel.empty()

    ch_fasta = file(params.fasta, checkIfExists: true)

    if (!params.sizes) {
        //
        // Create chromosome sizes file
        //
        CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
        ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    } else {
        ch_chrom_sizes = file(params.sizes)
    }

    ch_bam = Channel.fromFilePairs(params.bams, size: -1)
        .map {
            meta, fastq ->
            def fmeta = [:]
            // Set meta.id
            fmeta.id = meta
            // Set meta.single_end
            fmeta.single_end = true
            [ fmeta, fastq ]
        }

    SAMTOOLS_SORT ( ch_bam, ch_fasta )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    SAMTOOLS_SORT.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .join(SAMTOOLS_INDEX.out.csi, by: [0], remainder: true)
        .map {
            meta, bam, bai, csi ->
                if (bai) {
                    [ meta, bam, bai ]
                } else {
                    [ meta, bam, csi ]
                }
        }
        .set { ch_bam_bai }


    DREG_PREP (
        ch_bam_bai,
        ch_chrom_sizes,
        params.assay_type
    )

    DREG_RUN (
        DREG_PREP.out.plus_bigwig.join(DREG_PREP.out.minus_bigwig, by: [0]),
        params.dreg_model
    )
}
