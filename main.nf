#!/usr/bin/env nextflow

params.bams = null
params.sizes = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.dict"
params.fasta = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.assay_type = "groseq"
params.dreg_model = "https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata"

include { CUSTOM_GETCHROMSIZES } from './modules/nf-core/custom/getchromsizes/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'

include { DREG_PREP } from './modules/local/dreg_prep/main'
include { DREG_RUN } from './modules/local/dreg/main'

workflow {
    ch_chrom_sizes = Channel.empty()

    if (!params.sizes) {
        //
        // Create chromosome sizes file
        //
        ch_fasta = Channel.fromPath(params.fasta)
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
            if (fastq.size() == 1) {
                fmeta.single_end = true
            } else {
                fmeta.single_end = false
            }
            [ fmeta, fastq ]
        }

    SAMTOOLS_INDEX (
        ch_bam
    )

    ch_bam
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
