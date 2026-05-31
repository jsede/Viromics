process FASTQ_QC {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/qc" }, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    val min_qual

    output:
    tuple val(sample_id), path("${sample_id}_trimmed*.fastq.gz"), emit: trimmed_reads
    path "${sample_id}_fastp.json", emit: json_report
    path "${sample_id}_fastp.html", emit: html_report

    script:
    if (params.platform == 'illumina') {
        // Paired-end Illumina QC
        """
        # Trim adapters with BBDuk
        bbduk.sh \\
            in1=${reads[0]} in2=${reads[1]} \\
            out1=adapter_trimmed_R1.fastq.gz out2=adapter_trimmed_R2.fastq.gz \\
            ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo \\
            threads=${task.cpus}

        # Quality trim with fastp
        fastp \\
            -i adapter_trimmed_R1.fastq.gz -I adapter_trimmed_R2.fastq.gz \\
            -o ${sample_id}_trimmed_R1.fastq.gz -O ${sample_id}_trimmed_R2.fastq.gz \\
            -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \\
            --thread ${task.cpus} \\
            -q ${min_qual} \\
            --disable_adapter_trimming
        """
    } else {
        // Single-end ONT QC
        """
        # Trim adapters with BBDuk
        bbduk.sh \\
            in=${reads} \\
            out=adapter_trimmed.fastq.gz \\
            literal=AACGCTTTCGCGTTTTTCGTGCGCCGCTTCA ktrim=l k=18 mink=11 hdist=1 \\
            threads=${task.cpus}

        # Quality trim with fastp
        fastp \\
            -i adapter_trimmed.fastq.gz \\
            -o ${sample_id}_trimmed.fastq.gz \\
            -j ${sample_id}_fastp.json -h ${sample_id}_fastp.html \\
            --thread ${task.cpus} \\
            -q ${min_qual} \\
            --disable_adapter_trimming
        """
    }
}
