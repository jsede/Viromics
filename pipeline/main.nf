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
      --min_coverage        Minimum HQ depth for variant calling (default: 10 for Illumina, 50 for ONT)
      --min_freq            Minimum alt allele frequency to report variants (default: 0.2 for Illumina, 0.5 for ONT)
      --min_consensus_depth Minimum HQ depth to call consensus base (falls back to min_coverage)
      --min_consensus_freq  Minimum alt allele frequency for haploid consensus & major indels (default: 0.50)
      --min_map_qual        Minimum mapping quality (MAPQ) for alignment/depth/variants (default: 20)
      --min_base_qual       Minimum base quality for depth/variants (default: 20 for Illumina, 15 for ONT)
      --min_variant_qual    Minimum VCF QUAL score for variant calling (default: 20)
      --primer_tolerance    samtools ampliconclip tolerance in bp (default: 5 for Illumina, 15 for ONT)
      --min_clipped_len     Minimum read length retained after QC/clipping (default: 30 Illumina, 50 ONT)
      --ont_min_mean_qual   Minimum mean read Phred quality for ONT fastplong QC (default: 10)
      --ont_min_length      Minimum ONT read length for fastplong QC (default: 100)
      --ont_max_length      Maximum ONT read length for fastplong QC, 0 = unlimited (default: 0)
      --ont_adapter_fasta   Optional FASTA file of custom ONT adapter sequences (default: auto-detect)
      --disable_baq         Pass -B (--no-BAQ) to bcftools mpileup (default: true)
      --max_mpileup_depth   Max pileup depth cap; deeper contigs undergo deterministic strand-balanced downsampling (default: 100000)
      --downsample_seed     Fixed RNG seed for strand-balanced downsampling (default: 42)
      --use_iupac           Output IUPAC ambiguity codes for mixed biallelic SNPs in consensus (default: false)
      --min_iupac_freq      Lower AF bound (inclusive) for IUPAC SNP coding (default: 0.25)
      --max_iupac_freq      Upper AF bound (exclusive) for IUPAC SNP coding (default: 0.75)
      --outdir              Directory to write results to (default: './results')
    """
}

/*
========================================================================================
    DSL2 MODULE IMPORTS
========================================================================================
*/

include { FASTQ_QC } from './modules/fastq_qc.nf'
include { INDEX_REFERENCE; ALIGN_READS } from './modules/align_reads.nf'
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

    if (!(params.platform in ['illumina', 'ont'])) {
        error "Invalid --platform '${params.platform}'. Supported values are 'illumina' or 'ont'."
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

    def min_map_qual = (params.min_map_qual != null) ? params.min_map_qual : 20
    def min_base_qual = params.min_base_qual
    if (min_base_qual == null || min_base_qual == 'null' || min_base_qual == '') {
        min_base_qual = (params.min_qual != null && params.min_qual != 'null' && params.min_qual != '')
            ? params.min_qual
            : ((params.platform == 'ont') ? 15 : 20)
    }

    def primer_tolerance = params.primer_tolerance
    if (primer_tolerance == null || primer_tolerance == 'null' || primer_tolerance == '') {
        primer_tolerance = (params.platform == 'ont') ? 15 : 5
    }

    def min_clipped_len = params.min_clipped_len
    if (min_clipped_len == null || min_clipped_len == 'null' || min_clipped_len == '') {
        min_clipped_len = (params.platform == 'ont') ? 50 : 30
    }

    def min_variant_qual = (params.min_variant_qual != null) ? params.min_variant_qual : 20
    def min_consensus_freq = (params.min_consensus_freq != null) ? params.min_consensus_freq : 0.50
    def min_iupac_freq = (params.min_iupac_freq != null) ? params.min_iupac_freq : 0.25
    def max_iupac_freq = (params.max_iupac_freq != null) ? params.max_iupac_freq : 0.75
    def disable_baq = (params.disable_baq != null) ? params.disable_baq : true
    def max_mpileup_depth = (params.max_mpileup_depth != null) ? params.max_mpileup_depth : 100000
    def downsample_seed = (params.downsample_seed != null) ? params.downsample_seed : 42

    if (min_iupac_freq < 0 || min_iupac_freq >= max_iupac_freq || max_iupac_freq > 1) {
        error "Invalid IUPAC frequency bounds: require 0 <= min_iupac_freq (${min_iupac_freq}) < max_iupac_freq (${max_iupac_freq}) <= 1."
    }

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
    log.info "Min Consensus Freq  : ${min_consensus_freq}"
    log.info "Min Mapping Qual    : ${min_map_qual}"
    log.info "Min Base Qual       : ${min_base_qual}"
    log.info "Min Variant Qual    : ${min_variant_qual}"
    log.info "Primer Tolerance    : ${primer_tolerance} bp (min len: ${min_clipped_len} bp)"
    log.info "Disable BAQ (-B)    : ${disable_baq}"
    log.info "Max Mpileup Depth   : ${max_mpileup_depth} (seed: ${downsample_seed})"
    log.info "Use IUPAC           : ${params.use_iupac} (band: [${min_iupac_freq}, ${max_iupac_freq}))"
    if (params.platform == 'ont') {
        log.info "ONT Mean Qual / Len : >= Q${params.ont_min_mean_qual}, len >= ${params.ont_min_length} (max: ${params.ont_max_length ?: 'unlimited'})"
    }
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
    ont_adapter_fasta_ch = params.ont_adapter_fasta ? file(params.ont_adapter_fasta, checkIfExists: true) : []
    
    // Core Workflow Execution
    FASTQ_QC(
        reads_ch,
        params.platform,
        min_base_qual,
        min_clipped_len,
        params.ont_min_mean_qual,
        params.ont_min_length,
        params.ont_max_length,
        params.ont_adapter ?: '',
        ont_adapter_fasta_ch
    )
    if (params.platform == 'illumina') {
        INDEX_REFERENCE(ref_ch)
        bwa_index_ch = INDEX_REFERENCE.out.index.map { it[1] }
    } else {
        bwa_index_ch = Channel.value([])
    }
    ALIGN_READS(FASTQ_QC.out.trimmed_reads, ref_ch, bwa_index_ch, min_map_qual, params.platform)
    PRIMER_TRIM(ALIGN_READS.out.bam, primer_bed_ch, primer_tolerance, min_clipped_len, params.platform)
    
    VARIANT_CALL(
        PRIMER_TRIM.out.trimmed_bam,
        ref_ch,
        min_coverage,
        min_freq,
        min_map_qual,
        min_base_qual,
        min_variant_qual,
        disable_baq,
        max_mpileup_depth,
        downsample_seed,
        min_consensus_freq,
        params.use_iupac,
        min_iupac_freq,
        max_iupac_freq,
        params.platform
    )
    
    // Join BAMs and consensus-filtered VCFs on sample_id so subclonal alleles never leak into haploid consensus
    consensus_input_ch = PRIMER_TRIM.out.trimmed_bam.join(VARIANT_CALL.out.consensus_vcf)
    
    CONSENSUS_GEN(consensus_input_ch, ref_ch, min_consensus_depth, min_map_qual, min_base_qual, params.use_iupac)
    
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
