process VARIANT_CALL {
    tag "$sample_id"
    publishDir { "${params.outdir}/${sample_id}/variants" }, mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference
    val min_coverage
    val min_freq
    val min_map_qual
    val min_base_qual
    val min_variant_qual
    val disable_baq
    val max_mpileup_depth
    val downsample_seed
    val min_consensus_freq
    val use_iupac
    val min_iupac_freq
    val max_iupac_freq
    val platform

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: vcf
    tuple val(sample_id), path("${sample_id}_consensus.vcf"), emit: consensus_vcf

    script:
    def platform_mpileup_opts = (platform == 'ont') ? '-X ont-sup' : '-X illumina'
    def baq_opt = disable_baq ? '-B' : ''
    def idepth_cap = (max_mpileup_depth > 0) ? max_mpileup_depth : 1000000
    def iupac_flag = use_iupac ? '--use-iupac' : ''
    """
    # Clean blank lines / CRLF in staged reference for bcftools mpileup & norm
    awk '{ sub(/\\r\$/, ""); if (NF > 0) print }' ${reference} > clean_ref.fasta

    # Deterministic strand-balanced downsampling if depth exceeds max_mpileup_depth
    downsample_bam_strand_balanced.py \\
        --bam ${bam} \\
        --out ${sample_id}_pileup.bam \\
        --max-depth ${max_mpileup_depth} \\
        --seed ${downsample_seed} \\
        --map-qual ${min_map_qual} \\
        --base-qual ${min_base_qual} \\
        --threads ${task.cpus}

    # Call variants with diploid likelihoods (captures Illumina subclonal alleles cleanly),
    # disable BAQ (-B) by default, emit FORMAT/AD & FORMAT/DP, and normalize/split multiallelics
    bcftools mpileup -Ou -f clean_ref.fasta ${platform_mpileup_opts} ${baq_opt} \\
        -d ${max_mpileup_depth} -L ${idepth_cap} \\
        -q ${min_map_qual} -Q ${min_base_qual} \\
        --ns UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \\
        -a FORMAT/AD,FORMAT/DP,INFO/AD,INFO/ADF,INFO/ADR \\
        ${sample_id}_pileup.bam | \\
        bcftools call -mv --ploidy 2 -Ou | \\
        bcftools norm -m -any -f clean_ref.fasta -Ob -o ${sample_id}_raw.bcf

    # Filter on QUAL, high-quality FORMAT/DP, and per-allele FORMAT/AD[0:1]/FORMAT/DP
    bcftools filter -i "QUAL>=${min_variant_qual} && FORMAT/DP>=${min_coverage} && FORMAT/DP>0 && (FORMAT/AD[0:1]/FORMAT/DP)>=${min_freq}" -Ov -o ${sample_id}_filtered.vcf ${sample_id}_raw.bcf

    # Decouple reported variants (>= min_freq) from consensus VCF (>= min_consensus_freq or IUPAC SNP band)
    prepare_consensus_vcf.py \\
        --input ${sample_id}_filtered.vcf \\
        --reported-out ${sample_id}.vcf \\
        --consensus-out ${sample_id}_consensus.vcf \\
        --min-freq ${min_freq} \\
        --min-consensus-freq ${min_consensus_freq} \\
        --min-coverage ${min_coverage} \\
        --min-variant-qual ${min_variant_qual} \\
        ${iupac_flag} \\
        --min-iupac-freq ${min_iupac_freq} \\
        --max-iupac-freq ${max_iupac_freq}
    """
}
