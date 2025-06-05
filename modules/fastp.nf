process FASTP {
    tag "$prefix"
    label "fastp"
    publishDir "results/fastp", mode: 'copy'

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
        --thread $task.cpus \
        -h ${prefix}_fastp.html \
        -j ${prefix}_fastp.json \
        --detect_adapter_for_pe
    """
}