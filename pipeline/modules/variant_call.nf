process VARIANT_CALL {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/variants" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    val min_coverage
    val min_freq
    val map_qual
    val base_qual
    val min_variant_qual

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf

    script:
    def ploidy_arg = params.use_iupac ? "" : "--ploidy 1"
    """
    # Call variants using bcftools mpileup and call
    # Set -d 8000 and -L 8000 to prevent silent dropping of indels at depth > 250 in high-depth amplicons
    bcftools mpileup -Ou -f ${reference} -d 8000 -L 8000 -q ${map_qual} -Q ${base_qual} ${bam} | \\
        bcftools call ${ploidy_arg} -mv -Ob -o ${sample_id}_raw.bcf
        
    # Filter variants by quality, depth, and allele frequency
    bcftools filter -i "QUAL>=${min_variant_qual} && DP>=${min_coverage} && (DP4[2]+DP4[3])/(DP4[0]+DP4[1]+DP4[2]+DP4[3])>=${min_freq}" -Ov -o ${sample_id}.vcf ${sample_id}_raw.bcf
    """
}
