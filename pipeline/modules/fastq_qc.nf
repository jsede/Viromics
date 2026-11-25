process FASTQ_QC {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/qc" }, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val platform
    val min_base_qual
    val min_clipped_len
    val ont_min_mean_qual
    val ont_min_length
    val ont_max_length
    val ont_trim_adapters
    val ont_adapter
    path ont_adapter_fasta

    output:
    tuple val(sample_id), path("${sample_id}_trimmed*.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.json", emit: json_report
    path "${sample_id}_fastp.html", emit: html_report
    tuple val(sample_id), path("versions_fastq_qc.yml"), emit: versions

    script:
    if (platform == 'illumina') {
        // Paired-end Illumina QC directly via fastp (preserves true raw_reads counts)
        """
        fastp \\
            -i ${reads[0]} -I ${reads[1]} \\
            -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz \\
            -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \\
            --thread ${task.cpus} \\
            --detect_adapter_for_pe \\
            --cut_tail \\
            --trim_poly_g \\
            -q ${min_base_qual} \\
            -l ${min_clipped_len}
        fastp --version 2>&1 | head -n 1 | sed 's/^fastp /fastp: /' > versions_fastq_qc.yml
        """
    } else {
        // Single-end ONT QC via fastplong (mean Phred quality + length limits; -A prevents false auto-trimming of conserved amplicon ends)
        def max_len_opt = (ont_max_length && ont_max_length > 0) ? "--length_limit ${ont_max_length}" : ""
        def adapter_disable_opt = ont_trim_adapters ? "" : "-A"
        if (ont_adapter_fasta) {
            """
            fastplong \\
                -i ${reads} \\
                -o ${sample_id}_trimmed.fastq.gz \\
                -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \\
                -w ${task.cpus} \\
                -m ${ont_min_mean_qual} \\
                -u 100 \\
                -l ${ont_min_length} \\
                ${max_len_opt} \\
                -a ${ont_adapter_fasta}
            fastplong --version 2>&1 | head -n 1 | sed 's/^fastplong /fastplong: /' > versions_fastq_qc.yml
            """
        } else if (ont_adapter) {
            """
            printf ">ont_adapter\\n%s\\n" "${ont_adapter}" > ont_adapter_tmp.fasta
            fastplong \\
                -i ${reads} \\
                -o ${sample_id}_trimmed.fastq.gz \\
                -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \\
                -w ${task.cpus} \\
                -m ${ont_min_mean_qual} \\
                -u 100 \\
                -l ${ont_min_length} \\
                ${max_len_opt} \\
                -a ont_adapter_tmp.fasta
            fastplong --version 2>&1 | head -n 1 | sed 's/^fastplong /fastplong: /' > versions_fastq_qc.yml
            """
        } else {
            """
            fastplong \\
                -i ${reads} \\
                -o ${sample_id}_trimmed.fastq.gz \\
                -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \\
                -w ${task.cpus} \\
                -m ${ont_min_mean_qual} \\
                -u 100 \\
                -l ${ont_min_length} \\
                ${max_len_opt} \\
                ${adapter_disable_opt}
            fastplong --version 2>&1 | head -n 1 | sed 's/^fastplong /fastplong: /' > versions_fastq_qc.yml
            """
        }
    }
}
