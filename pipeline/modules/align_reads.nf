process INDEX_REFERENCE {
    tag "$reference"

    input:
    path reference

    output:
    tuple path(reference), path("${reference}.*"), emit: index

    script:
    """
    bwa index ${reference}
    """
}

process ALIGN_READS {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/alignment" }, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path reference
    path bwa_index
    val min_map_qual
    val platform

    output:
    tuple val(sample_id), path("${sample_id}_raw.bam"), emit: bam

    script:
    if (platform == 'illumina') {
        """
        # Align paired-end reads (pre-computed index staged), drop unmapped/secondary/supplementary & low MAPQ
        bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} | \\
            samtools view -@ ${task.cpus} -b -F 0x904 -q ${min_map_qual} - | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}_raw.bam -
        """
    } else {
        """
        # Align single-end long reads without secondary alignments, drop unmapped/secondary/supplementary & low MAPQ
        minimap2 -ax map-ont --secondary=no -t ${task.cpus} ${reference} ${reads} | \\
            samtools view -@ ${task.cpus} -b -F 0x904 -q ${min_map_qual} - | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}_raw.bam -
        """
    }
}
