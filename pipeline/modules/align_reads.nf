process ALIGN_READS {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/alignment" }, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_raw.bam"), emit: bam

    script:
    if (params.platform == 'illumina') {
        """
        # Index reference
        bwa index ${reference}
        
        # Align paired-end reads
        bwa mem -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} | \\
            samtools view -@ ${task.cpus} -bS - | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}_raw.bam -
        """
    } else {
        """
        # Align single-end long reads
        minimap2 -ax map-ont -t ${task.cpus} ${reference} ${reads} | \\
            samtools view -@ ${task.cpus} -bS - | \\
            samtools sort -@ ${task.cpus} -o ${sample_id}_raw.bam -
        """
    }
}
