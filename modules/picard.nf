process PICARD_COLLECT_MULTIPLE_METRICS {
    tag "$prefix"
    label "picard"
    
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
    tag "$prefix"
    label "picard"
    
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