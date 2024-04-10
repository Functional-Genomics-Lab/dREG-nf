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
        .map {
            validateInputSamplesheet(it)
        }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs.flatten() ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs.flatten() ]
        }

params.bwa_index = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/BWAIndex/version0.6.0/"
params.chromInfo = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/chromInfo.txt.gz"
params.assay_type = "GROseq"
params.sizes = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.dict"
params.fasta = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
params.dreg_model = "https://dreg.dnasequence.org/themes/dreg/assets/file/asvm.gdm.6.6M.20170828.rdata"
params.outdir = "results"

include { CAT_FASTQ } from './modules/nf-core/cat/fastq/main'
include { GUNZIP } from './modules/nf-core/gunzip/main'
include { PROSEQ2 } from './modules/local/proseq2/main'

workflow {
    CAT_FASTQ ( ch_input.multiple ).reads
        .mix(ch_input.single)
        .set { ch_cat_fastq }

    ch_chromInfo = Channel.empty()
    if (params.chromInfo.endsWith('.gz')) {
        ch_chromInfo = GUNZIP ( [ [:], file(params.chromInfo) ] ).gunzip.map { it[1] }
    } else {
        ch_chromInfo = Channel.fromPath(file(params.chromInfo))
    }

    PROSEQ2 (
        ch_cat_fastq,
        params.bwa_index,
        ch_chromInfo,
        params.assay_type,
    )
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}