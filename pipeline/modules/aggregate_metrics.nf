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
    val min_coverage
    val min_freq

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
    """
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
        --min-coverage ${min_coverage} \\
        --min-freq ${min_freq} \\
        --output viromics_run_report.json
    """
}
