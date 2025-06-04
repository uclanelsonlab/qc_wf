process MULTIQC {
    tag "multiqc"
        
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