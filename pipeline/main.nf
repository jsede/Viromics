#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
========================================================================================
    Viromics: Amplicon WGS Consensus calling & QC pipeline
========================================================================================
*/

// Help message
def helpMessage() {
    log.info"""
    Viromics Pipeline - Usage:
    
    nextflow run pipeline/main.nf --platform [illumina|ont] [reads/samplesheet options] --reference [fasta] [options]
    
    Required Arguments (One of reads or sample_sheet is required):
      --platform            Sequencing platform: 'illumina' or 'ont'
      --reads               Path pattern to FastQ files (e.g. "data/*_R{1,2}.fastq.gz" or "data/*.fastq.gz")
      --sample_sheet        Path to CSV samplesheet (alternative to --reads, contains: sample_id,platform,fastq_1,fastq_2)
      --reference           Path to reference genome (.fasta)
      
    Optional Arguments:
      --gff                 Path to gene annotation file (.gff3) for variant effect annotation
      --primer_bed          Path to primer scheme BED file for amplicon trimming
      --amplicons_bed       Path to explicit amplicon coordinates BED file
      --metadata            Path to sample metadata CSV file
      --min_coverage        Minimum depth for variant calling (default: 10 for Illumina, 50 for ONT)
      --min_freq            Minimum alt allele frequency for variants (default: 0.2 for Illumina, 0.5 for ONT)
      --min_consensus_depth Minimum depth to call consensus base (falls back to min_coverage)
      --use_iupac           Output IUPAC ambiguity codes for mixed sites in consensus (default: false)
      --outdir              Directory to write results to (default: './results')
    """
}

/*
========================================================================================
    DSL2 MODULE IMPORTS
========================================================================================
*/

include { FASTQ_QC } from './modules/fastq_qc.nf'
include { ALIGN_READS } from './modules/align_reads.nf'
include { PRIMER_TRIM } from './modules/primer_trim.nf'
include { VARIANT_CALL } from './modules/variant_call.nf'
include { CONSENSUS_GEN } from './modules/consensus_gen.nf'
include { AGGREGATE_METRICS } from './modules/aggregate_metrics.nf'

/*
=======================================================================================
    WORKFLOW
=======================================================================================
*/

workflow {
    // Show help message if requested or if required params are missing (one of reads or sample_sheet is required)
    if (params.help || !params.platform || (!params.reads && !params.sample_sheet) || !params.reference) {
        helpMessage()
        exit 0
    }

    // Set platform-specific default parameters dynamically
    def min_coverage = params.min_coverage
    if (min_coverage == null || min_coverage == 'null' || min_coverage == '') {
        min_coverage = (params.platform == 'ont') ? 50 : 10
    }
    
    def min_freq = params.min_freq
    if (min_freq == null || min_freq == 'null' || min_freq == '') {
        min_freq = (params.platform == 'ont') ? 0.5 : 0.2
    }
    
    def min_consensus_depth = params.min_consensus_depth
    if (min_consensus_depth == null || min_consensus_depth == 'null' || min_consensus_depth == '') {
        min_consensus_depth = min_coverage
    }

    def min_qual = params.min_qual ?: 20

    log.info "Starting Viromics Pipeline..."
    log.info "Platform            : ${params.platform}"
    log.info "Reads               : ${params.reads ?: 'None'}"
    log.info "Sample Sheet CSV    : ${params.sample_sheet ?: 'None'}"
    log.info "Reference           : ${params.reference}"
    log.info "GFF (GenBank)       : ${params.gff ?: 'None'}"
    log.info "Primer BED          : ${params.primer_bed ?: 'None'}"
    log.info "Amplicons BED       : ${params.amplicons_bed ?: 'None'}"
    log.info "Metadata CSV        : ${params.metadata ?: 'None'}"
    log.info "Min Coverage (DP)   : ${min_coverage}"
    log.info "Min Allele Freq (AF): ${min_freq}"
    log.info "Min Consensus Depth : ${min_consensus_depth}"
    log.info "Use IUPAC           : ${params.use_iupac}"
    log.info "OutDir              : ${params.outdir}"
    
    // Define input channels
    if (params.sample_sheet) {
        samplesheet_file = file(params.sample_sheet, checkIfExists: true)
        
        if (params.platform == 'illumina') {
            reads_ch = Channel
                .fromPath(samplesheet_file)
                .splitCsv(header: true)
                .map { row ->
                    if (!row.sample_id || !row.fastq_1 || !row.fastq_2) {
                        error "Samplesheet Error: row must contain 'sample_id', 'fastq_1', and 'fastq_2' for Illumina paired-end runs. Row: ${row}"
                    }
                    def fq1 = file(row.fastq_1, checkIfExists: true)
                    def fq2 = file(row.fastq_2, checkIfExists: true)
                    return tuple(row.sample_id, [fq1, fq2])
                }
        } else {
            reads_ch = Channel
                .fromPath(samplesheet_file)
                .splitCsv(header: true)
                .map { row ->
                    if (!row.sample_id || !row.fastq_1) {
                        error "Samplesheet Error: row must contain 'sample_id' and 'fastq_1' for ONT runs. Row: ${row}"
                    }
                    def fq1 = file(row.fastq_1, checkIfExists: true)
                    return tuple(row.sample_id, fq1)
                }
        }
    } else {
        if (params.platform == 'illumina') {
            reads_ch = Channel
                .fromFilePairs(params.reads, checkIfExists: true)
        } else {
            reads_ch = Channel
                .fromPath(params.reads, checkIfExists: true)
                .map { file -> tuple(file.simpleName, file) }
        }
    }
    
    ref_ch = file(params.reference, checkIfExists: true)
    gff_ch = params.gff ? file(params.gff, checkIfExists: true) : []
    primer_bed_ch = params.primer_bed ? file(params.primer_bed, checkIfExists: true) : []
    amplicons_bed_ch = params.amplicons_bed ? file(params.amplicons_bed, checkIfExists: true) : []
    metadata_ch = params.metadata ? file(params.metadata, checkIfExists: true) : []
    
    // Core Workflow Execution
    FASTQ_QC(reads_ch, min_qual)
    ALIGN_READS(FASTQ_QC.out.trimmed_reads, ref_ch)
    PRIMER_TRIM(ALIGN_READS.out.bam, primer_bed_ch)
    
    VARIANT_CALL(PRIMER_TRIM.out.trimmed_bam, ref_ch, min_coverage, min_qual)
    
    // Join the BAMs and VCFs channels on sample_id key to ensure they match
    consensus_input_ch = PRIMER_TRIM.out.trimmed_bam.join(VARIANT_CALL.out.vcf)
    
    CONSENSUS_GEN(consensus_input_ch, ref_ch, min_consensus_depth)
    
    // Aggregate all sample results
    run_id = new Date().format("yyyy-MM-dd_HH-mm-ss")
    
    // Gather outputs as lists to pass to aggregation process
    fastp_reports = FASTQ_QC.out.json_report.collect()
    trimmed_bams = PRIMER_TRIM.out.trimmed_bam.flatMap { tuple -> [tuple[1], tuple[2]] }.collect()
    sample_vcfs = VARIANT_CALL.out.vcf.map { tuple -> tuple[1] }.collect()
    sample_consensuses = CONSENSUS_GEN.out.consensus.map { tuple -> tuple[1] }.collect()
    
    AGGREGATE_METRICS(
        run_id,
        fastp_reports,
        trimmed_bams,
        sample_vcfs,
        sample_consensuses,
        ref_ch,
        gff_ch,
        primer_bed_ch,
        amplicons_bed_ch,
        metadata_ch,
        min_coverage,
        min_freq
    )
}
