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
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, meta + [ single_end:true ], [ fastq_1 ] ]
                } else {
                    return [ meta.id, meta + [ single_end:false ], [ fastq_1, fastq_2 ] ]
                }
        }
        .groupTuple()

params.bwa_index = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/version0.6.0/"
params.chromInfo = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz"
params.assay_type = "GROseq"
params.sizes = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.dict"
params.fasta = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.dreg_model = "https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata"
params.outdir = "results"

include { PROSEQ2 } from './modules/local/proseq2/main'

workflow {
    PROSEQ2 (
        ch_input,
        params.bwa_index,
        params.chromInfo,
        params.assay_type,
    )
}