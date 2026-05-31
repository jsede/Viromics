process VARIANT_CALL {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/variants" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    val min_coverage
    val min_qual

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    """
    # Call variants using bcftools mpileup and call
    bcftools mpileup -Ou -f ${reference} -d 8000 -q ${min_qual} -Q ${min_qual} ${bam} | \\
        bcftools call -mv -Ob -o ${sample_id}_raw.bcf
        
    # Filter variants by quality and depth
    bcftools filter -i "QUAL>=${min_qual} && DP>=${min_coverage}" -Ov -o ${sample_id}.vcf ${sample_id}_raw.bcf
    """
}
