process CHECK_AND_PROCESS_ALIGNMENT {
    tag "$prefix"
    label "samtools_process"
    
    input:
    path aligned_file
    path fasta
    val prefix
    
    output:
    path "${prefix}.bam", emit: bam
    path "${prefix}.bam.bai", emit: bam_index
    
    script:
    def file_ext = aligned_file.toString().toLowerCase()
    if (file_ext.endsWith('.cram')) {
        """
        samtools view -@ $task.cpus -b -T ${fasta} ${aligned_file} > ${prefix}.bam
        samtools index -@ $task.cpus ${prefix}.bam
        """
    } else if (file_ext.endsWith('.bam')) {
        """
        cp ${aligned_file} ${prefix}.bam
        samtools index -@ $task.cpus ${prefix}.bam
        """
    } else {
        error "Unsupported file format. Please provide either a BAM or CRAM file."
    }
}

// Process to convert BAM to CRAM
process BAM_TO_CRAM {
    label 'samtools_cram'
    publishDir "results/cram", mode: 'copy'
    
    input:
    path bam
    path fasta
    val prefix
    
    output:
    path "${prefix}.cram", emit: cram
    path "${prefix}.cram.crai", emit: cram_index
    
    when:
    params.cram
    
    script:
    """
    samtools view -C -T ${fasta} ${bam} > ${prefix}.cram
    samtools index ${prefix}.cram
    """
}