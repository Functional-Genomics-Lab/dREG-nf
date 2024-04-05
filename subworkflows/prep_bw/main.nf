// TODO
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


    // DREG_PREP (
    //     ch_bam_bai,
    //     ch_chrom_sizes,
    //     params.assay_type
    // )
}
