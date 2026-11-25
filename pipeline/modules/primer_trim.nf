process PRIMER_TRIM {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/alignment" }, mode: 'copy'

    input:
    tuple val(sample_id), path(raw_bam)
    path primer_bed
    val primer_tolerance
    val min_clipped_len
    val platform

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.bam"), path("${sample_id}_trimmed.bam.bai"), emit: trimmed_bam
    tuple val(sample_id), path("versions_primer_trim.yml"), emit: versions

    script:
    if (primer_bed) {
        if (platform == 'illumina') {
            """
            # Strand-aware primer clipping on both ends with tolerance & post-clip mate-fixing for paired-end Illumina
            samtools ampliconclip --both-ends --strand --soft-clip --tolerance ${primer_tolerance} --filter-len ${min_clipped_len} --no-excluded -@ ${task.cpus} -b ${primer_bed} ${raw_bam} | \\
                samtools sort -n -@ ${task.cpus} - | \\
                samtools fixmate -@ ${task.cpus} - - | \\
                samtools sort -@ ${task.cpus} -o ${sample_id}_trimmed.bam -
            samtools index ${sample_id}_trimmed.bam
            samtools --version 2>&1 | head -n 1 | sed 's/^samtools /samtools: /' > versions_primer_trim.yml
            """
        } else {
            """
            # Strand-aware primer clipping on both ends with ONT indel tolerance & minimum clipped insert length
            samtools ampliconclip --both-ends --strand --soft-clip --tolerance ${primer_tolerance} --filter-len ${min_clipped_len} --no-excluded -@ ${task.cpus} -b ${primer_bed} ${raw_bam} | \\
                samtools sort -@ ${task.cpus} -o ${sample_id}_trimmed.bam -
            samtools index ${sample_id}_trimmed.bam
            samtools --version 2>&1 | head -n 1 | sed 's/^samtools /samtools: /' > versions_primer_trim.yml
            """
        }
    } else {
        """
        # No primer scheme provided: preserve staged input symlink and index cleanly
        cp -L ${raw_bam} ${sample_id}_trimmed.bam
        samtools index ${sample_id}_trimmed.bam
        samtools --version 2>&1 | head -n 1 | sed 's/^samtools /samtools: /' > versions_primer_trim.yml
        """
    }
}
