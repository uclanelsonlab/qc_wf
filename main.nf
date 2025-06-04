#!/usr/bin/env nextflow

/*
 * WGS QC Workflow
 * Author: George Carvalho
 * Email: gcarvalhoneto@mednet.ucla.edu
 * Description: Pipeline for WGS QC stats check
 */

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp.nf'
include { CHECK_AND_PROCESS_ALIGNMENT } from './modules/samtools.nf'
include { PICARD_COLLECT_MULTIPLE_METRICS; PICARD_COLLECT_WGS_METRICS } from './modules/picard.nf'
include { MULTIQC } from './modules/multiqc.nf'


// Input parameters
params.fastq_r1 = null
params.fastq_r2 = null
params.prefix = null
params.fasta = null
params.aligned_file = null
params.java_mem = "-Xmx100g"

// Docker containers
params.fastp_container = "staphb/fastp:latest"
params.picard_container = "broadinstitute/picard:latest"
params.multiqc_container = "ewels/multiqc:latest"
params.samtools_container = "quay.io/biocontainers/samtools:1.21--h96c455f_1"

// Validate required parameters
if (!params.fastq_r1 || !params.fastq_r2 || !params.prefix || !params.fasta || !params.aligned_file) {
    error "Missing required parameters. Please provide fastq_r1, fastq_r2, prefix, fasta, and aligned_file"
}

// Input channels
Channel
    .fromPath(params.fastq_r1)
    .ifEmpty { error "Cannot find any reads matching: ${params.fastq_r1}" }
    .set { ch_fastq_r1 }

Channel
    .fromPath(params.fastq_r2)
    .ifEmpty { error "Cannot find any reads matching: ${params.fastq_r2}" }
    .set { ch_fastq_r2 }

Channel
    .fromPath(params.fasta)
    .ifEmpty { error "Cannot find reference genome: ${params.fasta}" }
    .set { ch_fasta }

Channel
    .fromPath(params.aligned_file)
    .ifEmpty { error "Cannot find aligned file: ${params.aligned_file}" }
    .set { ch_aligned_file }

// Main workflow
workflow {
    // Run Fastp
    FASTP(params.fastq_r1, params.fastq_r2, params.prefix)
    
    // Process aligned file (BAM or CRAM)
    CHECK_AND_PROCESS_ALIGNMENT(params.aligned_file, params.fasta, params.prefix)
    
    // Run Picard tools
    PICARD_COLLECT_MULTIPLE_METRICS(
        params.fasta,
        CHECK_AND_PROCESS_ALIGNMENT.out.bam,
        CHECK_AND_PROCESS_ALIGNMENT.out.bam_index,
        params.prefix,
        params.java_mem
    )
    PICARD_COLLECT_WGS_METRICS(
        params.fasta,
        CHECK_AND_PROCESS_ALIGNMENT.out.bam,
        CHECK_AND_PROCESS_ALIGNMENT.out.bam_index,
        params.prefix,
        params.java_mem
    )
    
    // Run MultiQC
    MULTIQC(
        FASTP.out.fastp_html,
        FASTP.out.fastp_json,
        PICARD_COLLECT_MULTIPLE_METRICS.out.metrics_files,
        PICARD_COLLECT_WGS_METRICS.out.wgs_metrics,
        params.prefix
    )
} 