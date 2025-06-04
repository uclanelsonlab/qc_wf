#!/usr/bin/env nextflow

/*
 * WGS QC Workflow
 * Author: George Carvalho
 * Email: gcarvalhoneto@mednet.ucla.edu
 * Description: Pipeline for WGS QC stats check
 */

// Enable DSL 2 syntax
nextflow.enable.dsl = 2

// Input parameters
params.fastq_r1 = null
params.fastq_r2 = null
params.prefix = null
params.fasta = null
params.bam_file = null
params.bam_index = null
params.java_mem = "-Xmx100g"

// Docker containers
params.fastp_container = "staphb/fastp:latest"
params.picard_container = "broadinstitute/picard:latest"
params.multiqc_container = "ewels/multiqc:latest"

// Validate required parameters
if (!params.fastq_r1 || !params.fastq_r2 || !params.prefix || !params.fasta || !params.bam_file || !params.bam_index) {
    error "Missing required parameters. Please provide fastq_r1, fastq_r2, prefix, fasta, bam_file, and bam_index"
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
    .fromPath(params.bam_file)
    .ifEmpty { error "Cannot find BAM file: ${params.bam_file}" }
    .set { ch_bam }

Channel
    .fromPath(params.bam_index)
    .ifEmpty { error "Cannot find BAM index file: ${params.bam_index}" }
    .set { ch_bam_index }

// Process definitions
process FASTP {
    container "${params.fastp_container}"
    
    input:
    path fastq_r1
    path fastq_r2
    val prefix
    
    output:
    path "${prefix}_fastp.html", emit: fastp_html
    path "${prefix}_fastp.json", emit: fastp_json
    
    script:
    """
    fastp -i ${fastq_r1} -I ${fastq_r2} \
        -h ${prefix}_fastp.html \
        -j ${prefix}_fastp.json \
        --detect_adapter_for_pe
    """
}

process PICARD_COLLECT_MULTIPLE_METRICS {
    container "${params.picard_container}"
    
    input:
    path fasta
    path bam
    path bam_index
    val prefix
    val java_mem
    
    output:
    path "${prefix}.*.metrics", emit: metrics_files
    path "${prefix}.*.pdf", emit: pdf_files
    
    script:
    """
    java ${java_mem} -jar /usr/picard/picard.jar CollectMultipleMetrics \
        -R ${fasta} \
        -I ${bam} \
        -O ${prefix} \
        --PROGRAM CollectAlignmentSummaryMetrics \
        --PROGRAM CollectInsertSizeMetrics \
        --PROGRAM QualityScoreDistribution \
        --PROGRAM MeanQualityByCycle \
        --PROGRAM CollectBaseDistributionByCycle
    """
}

process PICARD_COLLECT_WGS_METRICS {
    container "${params.picard_container}"
    
    input:
    path fasta
    path bam
    path bam_index
    val prefix
    val java_mem
    
    output:
    path "${prefix}.wgs_metrics", emit: wgs_metrics
    
    script:
    """
    ls ${bam} ${bam_index} ${fasta}
    java ${java_mem} -jar /usr/picard/picard.jar CollectWgsMetrics -R ${fasta} -I ${bam} -O ${prefix}.wgs_metrics
    """
}

process MULTIQC {
    container "${params.multiqc_container}"
    
    publishDir "results/multiqc", mode: 'copy'
    
    input:
    path fastp_html
    path fastp_json
    path metrics_files
    path wgs_metrics
    val prefix
    
    output:
    path "multiqc_report.html", emit: multiqc_html
    path "multiqc_data", emit: multiqc_data
    
    script:
    """
    multiqc . \
        --filename ${prefix}_multiqc_report.html \
        --outdir .
    """
}

// Main workflow
workflow {
    // Run Fastp
    FASTP(params.fastq_r1, params.fastq_r2, params.prefix)
    
    // Run Picard tools
    PICARD_COLLECT_MULTIPLE_METRICS(params.fasta, params.bam_file, params.bam_index, params.prefix, params.java_mem)
    PICARD_COLLECT_WGS_METRICS(params.fasta, params.bam_file, params.bam_index, params.prefix, params.java_mem)
    
    // Run MultiQC
    MULTIQC(
        FASTP.out.fastp_html,
        FASTP.out.fastp_json,
        PICARD_COLLECT_MULTIPLE_METRICS.out.metrics_files,
        PICARD_COLLECT_WGS_METRICS.out.wgs_metrics,
        params.prefix
    )
} 