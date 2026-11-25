process PRIMER_TRIM {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/alignment" }, mode: 'copy'

    input:
    tuple val(sample_id), path(raw_bam)
    path primer_bed

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.bam"), path("${sample_id}_trimmed.bam.bai"), emit: trimmed_bam

    script:
    if (primer_bed) {
        """
        # Remove primers using samtools ampliconclip clipping both ends
        samtools ampliconclip --both-ends -@ ${task.cpus} -b ${primer_bed} ${raw_bam} | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}_trimmed.bam -
        samtools index ${sample_id}_trimmed.bam
        """
    } else {
        """
        # No primer scheme provided, just sort and index the raw BAM
        mv ${raw_bam} ${sample_id}_trimmed.bam
        samtools index ${sample_id}_trimmed.bam
        """
    }
}
