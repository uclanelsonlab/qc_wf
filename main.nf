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
include { CHECK_AND_PROCESS_ALIGNMENT; BAM_TO_CRAM } from './modules/samtools.nf'
include { PICARD_COLLECT_MULTIPLE_METRICS; PICARD_COLLECT_WGS_METRICS } from './modules/picard.nf'
include { QUALIMAP_BAMQC } from './modules/qualimap.nf'
include { MULTIQC } from './modules/multiqc.nf'

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
        params.prefix
    )
    PICARD_COLLECT_WGS_METRICS(
        params.fasta,
        CHECK_AND_PROCESS_ALIGNMENT.out.bam,
        CHECK_AND_PROCESS_ALIGNMENT.out.bam_index,
        params.prefix
    )
    // Run Qualimap
    QUALIMAP_BAMQC(
        CHECK_AND_PROCESS_ALIGNMENT.out.bam,
        params.prefix
    )
    
    if (params.cram) {
        // Convert to CRAM if requested
        BAM_TO_CRAM(
            CHECK_AND_PROCESS_ALIGNMENT.out.bam,
            params.fasta,
            params.prefix
        )
    }
    
    // Run MultiQC
    MULTIQC(
        FASTP.out.fastp_html,
        FASTP.out.fastp_json,
        PICARD_COLLECT_MULTIPLE_METRICS.out.metrics_files,
        PICARD_COLLECT_WGS_METRICS.out.wgs_metrics,
        QUALIMAP_BAMQC.out.results,
        params.prefix
    )
} 