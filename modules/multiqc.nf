process MULTIQC {
    tag "$prefix"
    label "multiqc"
        
    publishDir "results/multiqc", mode: 'copy'
    
    input:
    path fastp_html
    path fastp_json
    path metrics_files
    path wgs_metrics
    path qualimap_results
    val prefix
    
    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional: true, emit: plots
    path "versions.yml"        , emit: versions
    
    script:
    """
    cp -r $qualimap_results .
    multiqc . \
        --filename ${prefix}_multiqc_report.html \
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}