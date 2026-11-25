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
      --run_id              Deterministic run ID for report provenance (default: viromics_<sessionId>)
      --gff                 Path to gene annotation file (.gff3) for variant effect annotation
      --primer_bed          Path to primer scheme BED file for amplicon trimming
      --amplicons_bed       Path to explicit amplicon coordinates BED file
      --metadata            Path to sample metadata CSV file
      --cds_models          Path to curated CDS frameshift JSON models (default: resources/annotation_models/cds_models.json)
      --min_coverage        Minimum HQ depth for variant calling (default: 10 for Illumina, 50 for ONT)
      --min_freq            Minimum alt allele frequency to report variants (default: 0.2 for Illumina, 0.5 for ONT)
      --min_consensus_depth Minimum HQ depth to call consensus base (falls back to min_coverage)
      --min_consensus_freq  Minimum alt allele frequency for haploid consensus & major indels (default: 0.50)
      --min_map_qual        Minimum mapping quality (MAPQ) for alignment/depth/variants (default: 20)
      --min_base_qual       Minimum base quality for depth/variants (default: 20 for Illumina, 15 for ONT)
      --min_variant_qual    Minimum VCF QUAL score for variant calling (default: 20)
      --primer_tolerance    samtools ampliconclip tolerance in bp (default: 5 for Illumina, 15 for ONT)
      --min_clipped_len     Minimum read length retained after QC/clipping (default: 30 Illumina, 50 ONT)
      --ont_mpileup_config  bcftools mpileup -X config for ONT (default: 'ont-sup'; pass 'none' to disable)
      --ont_min_mean_qual   Minimum mean read Phred quality for ONT fastplong QC (default: 10)
      --ont_min_length      Minimum ONT read length for fastplong QC (default: 100)
      --ont_max_length      Maximum ONT read length for fastplong QC, 0 = unlimited (default: 0)
      --ont_trim_adapters   Enable de novo fastplong auto-adapter trimming on ONT (default: false; passes -A to protect amplicon ends)
      --ont_adapter_fasta   Optional FASTA file of custom ONT adapter sequences (default: skipped unless ont_trim_adapters=true)
      --disable_baq         Pass -B (--no-BAQ) to bcftools mpileup (default: true)
      --max_mpileup_depth   Max pileup depth cap; deeper contigs undergo deterministic strand-balanced downsampling (default: 100000)
      --downsample_seed     Fixed RNG seed for strand-balanced downsampling (default: 42)
      --use_iupac           Output IUPAC ambiguity codes for mixed biallelic SNPs in consensus (default: false)
      --min_iupac_freq      Lower AF bound (inclusive) for IUPAC SNP coding (default: 0.25)
      --max_iupac_freq      Upper AF bound (exclusive) for IUPAC SNP coding (default: 0.75)
      --qc_pass_depth       QC PASS min mean depth on callable target (default: 100)
      --qc_pass_breadth     QC PASS min callable breadth % at min_coverage (default: 95.0)
      --qc_pass_max_n       QC PASS max internal N % (default: 5.0)
      --qc_warn_depth       QC WARN min mean depth on callable target (default: 25)
      --qc_warn_breadth     QC WARN min callable breadth % at min_coverage (default: 80.0)
      --qc_warn_max_n       QC WARN max internal N % (default: 10.0)
      --container_image     Container image tag or Singularity .sif path (default: 'viromics:2.0.0')
      --max_cpus            Maximum CPUs permitted per process (default: 16)
      --max_memory          Maximum memory permitted per process (default: '32.GB')
      --max_time            Maximum wall-clock time permitted per process (default: '24.h')
      --outdir              Directory to write results to (default: './results')
    """
}

def validateParameters(Map allParams, File schemaFile) {
    if (!schemaFile.exists()) {
        error "Parameter schema file not found: ${schemaFile}"
    }
    def schema = new groovy.json.JsonSlurper().parse(schemaFile)
    def schemaProps = [:]
    if (schema.'$defs' instanceof Map) {
        schema.'$defs'.each { _secName, secDef ->
            if (secDef?.properties instanceof Map) {
                schemaProps.putAll(secDef.properties)
            }
        }
    }
    if (schema.properties instanceof Map) {
        schemaProps.putAll(schema.properties)
    }

    allParams.each { key, val ->
        def kStr = key.toString()
        if (!schemaProps.containsKey(kStr)) {
            def exactPrefix = schemaProps.keySet().findAll { cand -> cand.toString().startsWith(kStr) }.sort()
            def subMatch = schemaProps.keySet().findAll { cand -> cand.toString().contains(kStr) || kStr.contains(cand.toString()) }.sort()
            def looseMatch = schemaProps.keySet().findAll { cand -> kStr.length() >= 4 && cand.toString().startsWith(kStr.take(4)) }.sort()
            def suggestion = exactPrefix ?: (subMatch ?: looseMatch)
            def hint = suggestion ? " Did you mean '--${suggestion[0]}'?" : ""
            error "Parameter validation error: Unknown parameter '--${kStr}'.${hint} Run with --help or inspect nextflow_schema.json."
        }

        if (val == null || val == 'null' || (val instanceof String && val.trim() == '')) {
            return
        }

        def prop = schemaProps[kStr]
        def expectedType = prop?.type
        if (expectedType == 'boolean') {
            boolean ok = (val instanceof Boolean) || (val instanceof String && val.trim().toLowerCase() in ['true', 'false'])
            if (!ok) {
                error "Parameter validation error: '--${kStr}' expects a boolean (true/false), got '${val}'."
            }
        } else if (expectedType == 'integer') {
            boolean ok = (val instanceof Integer || val instanceof Long) || (val instanceof String && (val.trim() ==~ /^-?\d+$/))
            if (!ok) {
                error "Parameter validation error: '--${kStr}' expects an integer value, got '${val}'."
            }
        } else if (expectedType == 'number') {
            boolean ok = ((val instanceof Number) && !(val instanceof Boolean)) || (val instanceof String && val.trim().isDouble())
            if (!ok) {
                error "Parameter validation error: '--${kStr}' expects a numeric value, got '${val}'."
            }
        } else if (expectedType == 'string') {
            if (val instanceof Boolean) {
                error "Parameter validation error: '--${kStr}' expects a string value, got boolean '${val}' (did you omit the argument value?)."
            }
        }

        if (prop?.'enum' instanceof List) {
            def vStr = val.toString().trim()
            if (!(vStr in prop.'enum'.collect { it.toString() })) {
                error "Parameter validation error: '--${kStr}' value '${vStr}' must be one of ${prop.'enum'}."
            }
        }

        if (expectedType in ['integer', 'number']) {
            def numVal = new BigDecimal(val.toString().trim())
            if (prop?.minimum != null && numVal < new BigDecimal(prop.minimum.toString())) {
                error "Parameter validation error: '--${kStr}' (${numVal}) must be >= ${prop.minimum}."
            }
            if (prop?.maximum != null && numVal > new BigDecimal(prop.maximum.toString())) {
                error "Parameter validation error: '--${kStr}' (${numVal}) must be <= ${prop.maximum}."
            }
        }
    }
}

def resolveSamplePath(String rawPath, String sampleId, String colName, samplesheetParent) {
    def direct = file(rawPath)
    if (direct.exists()) return direct
    if (!new File(rawPath).isAbsolute()) {
        def rel = file("${samplesheetParent}/${rawPath}")
        if (rel.exists()) return rel
    }
    error "Samplesheet Error: file not found for sample '${sampleId}' (${colName}): '${rawPath}'"
}

def validateSamplesheetRow(Map row, Set seenSampleIds, String expectedPlatform, String sheetPath) {
    def sid = row.sample_id?.toString()?.trim()
    if (!sid) {
        error "Samplesheet Error: missing or empty 'sample_id' in row: ${row}"
    }
    if (!seenSampleIds.add(sid)) {
        error "Samplesheet Error: duplicate sample_id '${sid}' in ${sheetPath}."
    }
    def rowPlatform = row.platform?.toString()?.trim()?.toLowerCase()
    if (rowPlatform) {
        if (!(rowPlatform in ['illumina', 'ont'])) {
            error "Samplesheet Error: invalid platform '${rowPlatform}' for sample '${sid}'. Supported values are 'illumina' or 'ont'."
        }
        if (rowPlatform != expectedPlatform) {
            error "Samplesheet Error: sample '${sid}' declares platform '${rowPlatform}' which conflicts with --platform '${expectedPlatform}'. Viromics processes one platform per run; split mixed-platform samplesheets or pass matching --platform."
        }
    }
    return sid
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
    // Validate CLI & config parameters against nextflow_schema.json before executing any tasks
    validateParameters(params, new File("${projectDir}/../nextflow_schema.json"))

    // Show help message if requested or if required params are missing (one of reads or sample_sheet is required)
    if (params.help || !params.platform || (!params.reads && !params.sample_sheet) || !params.reference) {
        helpMessage()
        exit 0
    }

    if (!(params.platform in ['illumina', 'ont'])) {
        error "Invalid --platform '${params.platform}'. Supported values are 'illumina' or 'ont'."
    }

    // Set platform-specific default parameters dynamically with strict numeric/boolean coercion
    def min_coverage = (params.min_coverage == null || params.min_coverage == 'null' || params.min_coverage == '')
        ? ((params.platform == 'ont') ? 50 : 10)
        : (params.min_coverage as int)

    def min_freq = (params.min_freq == null || params.min_freq == 'null' || params.min_freq == '')
        ? ((params.platform == 'ont') ? 0.5d : 0.2d)
        : (params.min_freq as double)

    def min_consensus_depth = (params.min_consensus_depth == null || params.min_consensus_depth == 'null' || params.min_consensus_depth == '')
        ? min_coverage
        : (params.min_consensus_depth as int)

    def min_map_qual = (params.min_map_qual != null) ? (params.min_map_qual as int) : 20
    def min_base_qual = (params.min_base_qual == null || params.min_base_qual == 'null' || params.min_base_qual == '')
        ? ((params.min_qual != null && params.min_qual != 'null' && params.min_qual != '')
            ? (params.min_qual as int)
            : ((params.platform == 'ont') ? 15 : 20))
        : (params.min_base_qual as int)

    def primer_tolerance = (params.primer_tolerance == null || params.primer_tolerance == 'null' || params.primer_tolerance == '')
        ? ((params.platform == 'ont') ? 15 : 5)
        : (params.primer_tolerance as int)

    def min_clipped_len = (params.min_clipped_len == null || params.min_clipped_len == 'null' || params.min_clipped_len == '')
        ? ((params.platform == 'ont') ? 50 : 30)
        : (params.min_clipped_len as int)

    def run_id = (params.run_id != null && params.run_id != 'null' && params.run_id.toString().trim() != '')
        ? params.run_id.toString().trim()
        : "viromics_${workflow.sessionId.toString().take(8)}"
    def ont_mpileup_config = (params.ont_mpileup_config != null && params.ont_mpileup_config != 'null')
        ? params.ont_mpileup_config.toString().trim()
        : 'ont-sup'
    def ont_min_mean_qual = (params.ont_min_mean_qual != null) ? (params.ont_min_mean_qual as int) : 10
    def ont_min_length = (params.ont_min_length != null) ? (params.ont_min_length as int) : 100
    def ont_max_length = (params.ont_max_length != null) ? (params.ont_max_length as int) : 0
    def ont_trim_adapters = (params.ont_trim_adapters in [true, 'true', 'True', 'TRUE'])

    def min_variant_qual = (params.min_variant_qual != null) ? (params.min_variant_qual as int) : 20
    def min_consensus_freq = (params.min_consensus_freq != null) ? (params.min_consensus_freq as double) : 0.50d
    def min_iupac_freq = (params.min_iupac_freq != null) ? (params.min_iupac_freq as double) : 0.25d
    def max_iupac_freq = (params.max_iupac_freq != null) ? (params.max_iupac_freq as double) : 0.75d
    def disable_baq = (params.disable_baq == null) ? true : (params.disable_baq in [true, 'true', 'True', 'TRUE'])
    def max_mpileup_depth = (params.max_mpileup_depth != null) ? (params.max_mpileup_depth as int) : 100000
    def downsample_seed = (params.downsample_seed != null) ? (params.downsample_seed as int) : 42
    def use_iupac = (params.use_iupac in [true, 'true', 'True', 'TRUE'])

    def qc_pass_depth = (params.qc_pass_depth != null) ? (params.qc_pass_depth as int) : 100
    def qc_pass_breadth = (params.qc_pass_breadth != null) ? (params.qc_pass_breadth as double) : 95.0d
    def qc_pass_max_n = (params.qc_pass_max_n != null) ? (params.qc_pass_max_n as double) : 5.0d
    def qc_warn_depth = (params.qc_warn_depth != null) ? (params.qc_warn_depth as int) : 25
    def qc_warn_breadth = (params.qc_warn_breadth != null) ? (params.qc_warn_breadth as double) : 80.0d
    def qc_warn_max_n = (params.qc_warn_max_n != null) ? (params.qc_warn_max_n as double) : 10.0d

    if (min_iupac_freq < 0 || min_iupac_freq >= max_iupac_freq || max_iupac_freq > 1) {
        error "Invalid IUPAC frequency bounds: require 0 <= min_iupac_freq (${min_iupac_freq}) < max_iupac_freq (${max_iupac_freq}) <= 1."
    }
    if (qc_warn_depth < 0 || qc_pass_depth < qc_warn_depth) {
        error "Invalid QC depth thresholds: require 0 <= qc_warn_depth (${qc_warn_depth}) <= qc_pass_depth (${qc_pass_depth})."
    }
    if (qc_warn_breadth < 0 || qc_warn_breadth > qc_pass_breadth || qc_pass_breadth > 100) {
        error "Invalid QC breadth thresholds: require 0 <= qc_warn_breadth (${qc_warn_breadth}) <= qc_pass_breadth (${qc_pass_breadth}) <= 100."
    }
    if (qc_pass_max_n < 0 || qc_pass_max_n > qc_warn_max_n || qc_warn_max_n > 100) {
        error "Invalid QC max N% thresholds: require 0 <= qc_pass_max_n (${qc_pass_max_n}) <= qc_warn_max_n (${qc_warn_max_n}) <= 100."
    }

    log.info "Starting Viromics Pipeline..."
    log.info "Run ID              : ${run_id}"
    log.info "Platform            : ${params.platform}"
    log.info "Reads               : ${params.reads ?: 'None'}"
    log.info "Sample Sheet CSV    : ${params.sample_sheet ?: 'None'}"
    log.info "Reference           : ${params.reference}"
    log.info "GFF (GenBank)       : ${params.gff ?: 'None'}"
    log.info "Primer BED          : ${params.primer_bed ?: 'None'}"
    log.info "Amplicons BED       : ${params.amplicons_bed ?: 'None'}"
    log.info "Metadata CSV        : ${params.metadata ?: 'None'}"
    log.info "Curated CDS Models  : ${params.cds_models ?: 'None'}"
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
    log.info "Use IUPAC           : ${use_iupac} (band: [${min_iupac_freq}, ${max_iupac_freq}))"
    log.info "QC Tiers (PASS/WARN): depth>=${qc_pass_depth}/${qc_warn_depth}, breadth>=${qc_pass_breadth}%/${qc_warn_breadth}%, maxN<=${qc_pass_max_n}%/${qc_warn_max_n}%"
    if (params.platform == 'ont') {
        log.info "ONT Mean Qual / Len : >= Q${ont_min_mean_qual}, len >= ${ont_min_length} (max: ${ont_max_length ?: 'unlimited'}), trim_adapters=${ont_trim_adapters}"
        log.info "ONT mpileup config  : ${ont_mpileup_config ?: 'none'}"
    }
    log.info "OutDir              : ${params.outdir}"
    
    // Define input channels
    if (params.sample_sheet) {
        samplesheet_file = file(params.sample_sheet, checkIfExists: true)
        def samplesheetParent = samplesheet_file.toAbsolutePath().parent
        def seenSampleIds = Collections.synchronizedSet([] as Set)

        if (params.platform == 'illumina') {
            reads_ch = Channel
                .fromPath(samplesheet_file)
                .splitCsv(header: true)
                .map { row ->
                    def sid = validateSamplesheetRow(row, seenSampleIds, params.platform, params.sample_sheet)
                    def r1 = row.fastq_1?.toString()?.trim()
                    def r2 = row.fastq_2?.toString()?.trim()
                    if (!r1 || !r2) {
                        error "Samplesheet Error: row for '${sid}' must contain 'fastq_1' and 'fastq_2' for Illumina paired-end runs. Row: ${row}"
                    }
                    def fq1 = resolveSamplePath(r1, sid, 'fastq_1', samplesheetParent)
                    def fq2 = resolveSamplePath(r2, sid, 'fastq_2', samplesheetParent)
                    return tuple(sid, [fq1, fq2])
                }
        } else {
            reads_ch = Channel
                .fromPath(samplesheet_file)
                .splitCsv(header: true)
                .map { row ->
                    def sid = validateSamplesheetRow(row, seenSampleIds, params.platform, params.sample_sheet)
                    def r1 = row.fastq_1?.toString()?.trim()
                    if (!r1) {
                        error "Samplesheet Error: row for '${sid}' must contain 'fastq_1' for ONT runs. Row: ${row}"
                    }
                    def fq1 = resolveSamplePath(r1, sid, 'fastq_1', samplesheetParent)
                    return tuple(sid, fq1)
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
    cds_models_ch = params.cds_models ? file(params.cds_models, checkIfExists: true) : []
    
    // Core Workflow Execution
    FASTQ_QC(
        reads_ch,
        params.platform,
        min_base_qual,
        min_clipped_len,
        ont_min_mean_qual,
        ont_min_length,
        ont_max_length,
        ont_trim_adapters,
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
        use_iupac,
        min_iupac_freq,
        max_iupac_freq,
        params.platform,
        ont_mpileup_config
    )
    
    // Join BAMs and consensus-filtered VCFs on sample_id so subclonal alleles never leak into haploid consensus
    consensus_input_ch = PRIMER_TRIM.out.trimmed_bam.join(VARIANT_CALL.out.consensus_vcf)
    
    CONSENSUS_GEN(consensus_input_ch, ref_ch, min_consensus_depth, min_map_qual, min_base_qual, use_iupac)
    
    // Select one deterministic versions_*.yml per process (lexicographically first sample_id) so -resume caches AGGREGATE_METRICS
    ch_versions = FASTQ_QC.out.versions.toSortedList { a, b -> a[0] <=> b[0] }.map { it[0][1] }
        .mix(ALIGN_READS.out.versions.toSortedList { a, b -> a[0] <=> b[0] }.map { it[0][1] })
        .mix(PRIMER_TRIM.out.versions.toSortedList { a, b -> a[0] <=> b[0] }.map { it[0][1] })
        .mix(VARIANT_CALL.out.versions.toSortedList { a, b -> a[0] <=> b[0] }.map { it[0][1] })
        .mix(CONSENSUS_GEN.out.versions.toSortedList { a, b -> a[0] <=> b[0] }.map { it[0][1] })
    if (params.platform == 'illumina') {
        ch_versions = ch_versions.mix(INDEX_REFERENCE.out.versions)
    }
    versions_files = ch_versions.toSortedList { a, b -> a.name <=> b.name }

    // Deterministic run provenance (no wall-clock timestamps so -resume caches AGGREGATE_METRICS)
    def run_config_map = [
        run_id: run_id,
        platform: params.platform,
        min_coverage: min_coverage,
        min_freq: min_freq,
        min_consensus_depth: min_consensus_depth,
        min_consensus_freq: min_consensus_freq,
        min_map_qual: min_map_qual,
        min_base_qual: min_base_qual,
        min_variant_qual: min_variant_qual,
        primer_tolerance: primer_tolerance,
        min_clipped_len: min_clipped_len,
        ont_mpileup_config: (params.platform == 'ont') ? (ont_mpileup_config ?: 'none') : null,
        ont_min_mean_qual: (params.platform == 'ont') ? ont_min_mean_qual : null,
        ont_min_length: (params.platform == 'ont') ? ont_min_length : null,
        ont_max_length: (params.platform == 'ont') ? ont_max_length : null,
        ont_trim_adapters: (params.platform == 'ont') ? ont_trim_adapters : null,
        disable_baq: disable_baq,
        max_mpileup_depth: max_mpileup_depth,
        downsample_seed: downsample_seed,
        use_iupac: use_iupac,
        min_iupac_freq: min_iupac_freq,
        max_iupac_freq: max_iupac_freq,
        reference: file(params.reference).name,
        gff: params.gff ? file(params.gff).name : null,
        primer_bed: params.primer_bed ? file(params.primer_bed).name : null,
        amplicons_bed: params.amplicons_bed ? file(params.amplicons_bed).name : null,
        metadata: params.metadata ? file(params.metadata).name : null,
        cds_models: params.cds_models ? file(params.cds_models).name : null,
        nextflow_version: workflow.nextflow.version.toString(),
        profile: workflow.profile.toString(),
        qc_thresholds: [
            pass_depth: qc_pass_depth,
            pass_breadth: qc_pass_breadth,
            pass_max_n: qc_pass_max_n,
            warn_depth: qc_warn_depth,
            warn_breadth: qc_warn_breadth,
            warn_max_n: qc_warn_max_n
        ]
    ]
    def run_config_json = groovy.json.JsonOutput.toJson(run_config_map)
    
    // Gather outputs in deterministic filename order so -resume reliably caches AGGREGATE_METRICS
    fastp_reports = FASTQ_QC.out.json_report.toSortedList { a, b -> a.name <=> b.name }
    trimmed_bams = PRIMER_TRIM.out.trimmed_bam.flatMap { tuple -> [tuple[1], tuple[2]] }.toSortedList { a, b -> a.name <=> b.name }
    sample_vcfs = VARIANT_CALL.out.vcf.map { tuple -> tuple[1] }.toSortedList { a, b -> a.name <=> b.name }
    sample_consensuses = CONSENSUS_GEN.out.consensus.map { tuple -> tuple[1] }.toSortedList { a, b -> a.name <=> b.name }
    
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
        cds_models_ch,
        versions_files,
        run_config_json,
        min_coverage,
        min_freq,
        min_consensus_depth,
        min_consensus_freq,
        min_map_qual,
        min_base_qual,
        qc_pass_depth,
        qc_pass_breadth,
        qc_pass_max_n,
        qc_warn_depth,
        qc_warn_breadth,
        qc_warn_max_n
    )
}
