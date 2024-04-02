#!/usr/bin/env nextflow

params.bams = null
params.sizes = null
params.assay_type = "groseq"
params.dreg_model = "https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata"

include { CUSTOM_GETCHROMSIZES } from './modules/nf-core/custom/getchromsizes/main'
include { DREG_PREP } from './modules/local/dreg_prep/main'
include { DREG_RUN } from './modules/local/dreg/main'

workflow {
    ch_chrom_sizes = Channel.empty()

    if (!params.sizes) {
        //
        // Create chromosome sizes file
        //
        CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
        ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
        ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    } else {
        ch_chrom_sizes = params.sizes
    }

    DREG_PREP (
        params.bams,
        ch_chrom_sizes,
        params.assay_type
    )

    DREG_RUN (
        DREG_PREP.out.plus_bigwig.join(DREG_PREP.out.minus_bigwig, by: [0]),
        params.dreg_model
    )
}