process CONSENSUS_GEN {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/consensus" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(vcf)
    path reference
    val min_consensus_depth

    output:
    tuple val(sample_id), path("${sample_id}_consensus.fasta"), emit: consensus

    script:
    def hap_arg = params.use_iupac ? "-H I" : "-H A"
    """
    # Clean empty lines from reference FASTA to prevent bcftools consensus from skipping contigs
    grep -v '^\$' ${reference} > clean_reference.fasta

    # Create coverage depth mask (regions below min_consensus_depth, including unmapped 0x coverage contigs)
    # -J includes reads with deletions so genuine biological deletions are not falsely masked with N
    samtools depth -aa -J ${bam} | \\
        awk '\$3 < ${min_consensus_depth} {print \$1 "\\t" \$2-1 "\\t" \$2}' > low_coverage_mask.bed
        
    # Compress and index VCF
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    
    # Generate consensus calling bcftools consensus, mask low-cov sites with N
    bcftools consensus ${hap_arg} -f clean_reference.fasta -m low_coverage_mask.bed ${vcf}.gz > ${sample_id}_consensus.fasta
    
    # Replace header to prefix with sample ID and retain the reference name in a portable way
    awk '/^>/ {print ">${sample_id}_" substr(\$0, 2); next} {print}' ${sample_id}_consensus.fasta > tmp_consensus.fasta && mv tmp_consensus.fasta ${sample_id}_consensus.fasta
    """
}
