process CHECK_AND_PROCESS_ALIGNMENT {
    tag "$prefix"
    label "samtools"
    
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
        samtools view -b -T ${fasta} ${aligned_file} > ${prefix}.bam
        samtools index ${prefix}.bam
        """
    } else if (file_ext.endsWith('.bam')) {
        """
        cp ${aligned_file} ${prefix}.bam
        samtools index ${prefix}.bam
        """
    } else {
        error "Unsupported file format. Please provide either a BAM or CRAM file."
    }
}