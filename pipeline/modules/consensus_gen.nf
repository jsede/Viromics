process CONSENSUS_GEN {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/consensus" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai), path(consensus_vcf)
    path reference
    val min_consensus_depth
    val min_map_qual
    val min_base_qual
    val use_iupac

    output:
    tuple val(sample_id), path("${sample_id}_consensus.fasta"), emit: consensus
    tuple val(sample_id), path("versions_consensus_gen.yml"), emit: versions

    script:
    def hap_arg = use_iupac ? "-H I" : "-H A"
    """
    # Clean CRLF and empty lines from reference FASTA to prevent bcftools consensus from skipping contigs
    awk '{ sub(/\\r\$/, ""); if (NF > 0) print }' ${reference} > clean_reference.fasta

    # Create coverage depth mask synced with VARIANT_CALL MAPQ/BaseQ/flag filters
    # -J includes reads with deletions so genuine biological deletions are not falsely masked with N
    samtools depth -aa -J -Q ${min_map_qual} -q ${min_base_qual} -G 0xF04 ${bam} | \\
        awk -v min_dp="${min_consensus_depth}" '\$3 < min_dp {print \$1 "\\t" \$2-1 "\\t" \$2}' > low_coverage_mask.bed
        
    # Compress and index consensus VCF
    bgzip -c ${consensus_vcf} > ${consensus_vcf}.gz
    tabix -p vcf ${consensus_vcf}.gz
    
    # Generate consensus calling bcftools consensus, mask low-cov sites with N
    bcftools consensus ${hap_arg} -f clean_reference.fasta -m low_coverage_mask.bed ${consensus_vcf}.gz > ${sample_id}_consensus.fasta
    
    # Replace header to prefix with sample ID and retain clean contig name (strip FASTA description)
    awk '/^>/ {split(substr(\$0, 2), a, /[ \\t]+/); print ">${sample_id}_" a[1]; next} {print}' ${sample_id}_consensus.fasta > tmp_consensus.fasta && mv tmp_consensus.fasta ${sample_id}_consensus.fasta

    samtools --version 2>&1 | head -n 1 | sed 's/^samtools /samtools: /' > versions_consensus_gen.yml
    bcftools --version 2>&1 | head -n 1 | sed 's/^bcftools /bcftools: /' >> versions_consensus_gen.yml
    (tabix --version 2>&1 || true) | head -n 1 | sed 's/^tabix[^0-9]*/tabix: /' >> versions_consensus_gen.yml
    """
}
