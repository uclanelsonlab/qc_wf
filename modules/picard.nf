process PICARD_COLLECT_MULTIPLE_METRICS {
    tag "$prefix"
    label "picard"
    
    input:
    path fasta
    path bam
    path bam_index
    val prefix
    
    output:
    path "${prefix}.*.metrics", emit: metrics_files
    path "${prefix}.*.pdf", emit: pdf_files
    
    script:
    def avail_mem = 50000
    if (!task.memory) {
        log.info '[Picard CollectMultipleMetrics] Available memory not known - defaulting to 50GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    java -Xmx${avail_mem}M -jar /usr/picard/picard.jar CollectMultipleMetrics \
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
    tag "$prefix"
    label "picard"
    
    input:
    path fasta
    path bam
    path bam_index
    val prefix
    
    output:
    path "${prefix}.wgs_metrics", emit: wgs_metrics
    
    script:
    def avail_mem = 50000
    if (!task.memory) {
        log.info '[Picard CollectWgsMetrics] Available memory not known - defaulting to 50GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    ls ${bam} ${bam_index} ${fasta}
    java -Xmx${avail_mem}M -jar /usr/picard/picard.jar CollectWgsMetrics -R ${fasta} -I ${bam} -O ${prefix}.wgs_metrics
    """
}