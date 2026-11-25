process AGGREGATE_METRICS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val run_id
    path fastp_jsons
    path bams
    path vcfs
    path consensuses
    path reference
    path gff
    path primer_bed
    path amplicons_bed
    path metadata
    path cds_models
    path versions_files
    val run_config_json
    val min_coverage
    val min_freq
    val min_consensus_depth
    val min_consensus_freq
    val min_map_qual
    val min_base_qual
    val qc_pass_depth
    val qc_pass_breadth
    val qc_pass_max_n
    val qc_warn_depth
    val qc_warn_breadth
    val qc_warn_max_n

    output:
    path "viromics_run_report.json"

    script:
    def fastp_list = fastp_jsons.join(',')
    def bam_list = bams.findAll { it.name.endsWith('.bam') }.join(',')
    def vcf_list = vcfs.join(',')
    def consensus_list = consensuses.join(',')
    def gff_arg = gff ? "--gff ${gff}" : ""
    def primers_arg = primer_bed ? "--primers ${primer_bed}" : ""
    def amplicons_bed_arg = amplicons_bed ? "--amplicons-bed ${amplicons_bed}" : ""
    def metadata_arg = metadata ? "--metadata ${metadata}" : ""
    def cds_models_arg = cds_models ? "--cds-models ${cds_models}" : ""
    def versions_list = (versions_files instanceof List ? versions_files : [versions_files]).findAll { it }.collect { it.toString() }.join(',')
    def versions_arg = versions_list ? "--versions \"${versions_list}\"" : ""
    """
    cat <<'EOF' > run_config.json
${run_config_json}
EOF

    parse_metrics.py \\
        --run-id "${run_id}" \\
        --fastp "${fastp_list}" \\
        --bams "${bam_list}" \\
        --vcfs "${vcf_list}" \\
        --consensuses "${consensus_list}" \\
        --reference "${reference}" \\
        ${gff_arg} \\
        ${primers_arg} \\
        ${amplicons_bed_arg} \\
        ${metadata_arg} \\
        ${cds_models_arg} \\
        ${versions_arg} \\
        --run-config run_config.json \\
        --min-coverage ${min_coverage} \\
        --min-freq ${min_freq} \\
        --min-consensus-depth ${min_consensus_depth} \\
        --min-consensus-freq ${min_consensus_freq} \\
        --min-map-qual ${min_map_qual} \\
        --min-base-qual ${min_base_qual} \\
        --qc-pass-depth ${qc_pass_depth} \\
        --qc-pass-breadth ${qc_pass_breadth} \\
        --qc-pass-max-n ${qc_pass_max_n} \\
        --qc-warn-depth ${qc_warn_depth} \\
        --qc-warn-breadth ${qc_warn_breadth} \\
        --qc-warn-max-n ${qc_warn_max_n} \\
        --output viromics_run_report.json
    """
}
