#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline
if (params.help) {
   log.info paramsHelp("nextflow run funlab/dREG-nf --input input_file.csv")
   exit 0
}

// TODO Validate input parameters
// validateParameters()

// Print summary of supplied parameters
log.info paramsSummaryLog(workflow)

// Create a new channel of metadata from a sample sheet
// NB: `input` corresponds to `params.input` and associated sample sheet schema
ch_input = Channel.fromSamplesheet("input")

params.sizes = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.dict"
params.fasta = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.assay_type = "GROseq"
params.dreg_model = "https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata"

workflow {
    PROSEQ2 (

    )
}