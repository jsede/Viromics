import fs from 'node:fs';
import path from 'node:path';
import { describe, expect, it } from 'vitest';
import {
  computeStableReportId,
  deduplicateConsequences,
  formatConsequences,
  formatDepth,
  formatMetric,
  getBestReference,
  getMappedReads,
  getSampleStatus,
  validateUploadJSON,
} from '../reportSchema.js';

describe('reportSchema v2 validation & v1 fallback', () => {
  it('validates v2 report with nullable ONT duplication_rate and null distances', () => {
    const v2Report = {
      schema_version: '2.0.0',
      run_id: 'ont_run_1',
      min_coverage_threshold: 50,
      run_config: {
        platform: 'ont',
        min_coverage: 50,
        min_freq: 0.5,
        qc_thresholds: {
          pass_depth: 100,
          pass_breadth: 95,
          pass_max_n: 5,
          warn_depth: 25,
          warn_breadth: 80,
          warn_max_n: 10,
        },
      },
      references: [{ name: 'FIPV_NC_002306', length: 434, sequence: 'A'.repeat(434) }],
      distance_matrices: {
        FIPV_NC_002306: [
          [0, null],
          [null, 0],
        ],
      },
      samples: [
        {
          sample_id: 'barcode01',
          qc_status: 'PASS',
          qc_reasons: ['meets_all_pass_thresholds'],
          primary_reference: 'FIPV_NC_002306',
          qc: {
            raw_reads: 1000,
            raw_bases: 400000,
            trimmed_reads: 980,
            trimmed_bases: 390000,
            q30_rate_raw: 0.52,
            q30_rate_trimmed: 0.56,
            duplication_rate: null,
          },
          coverage: {
            FIPV_NC_002306: {
              mean_depth: 800,
              callable_mean_depth: 880,
              coverage_breadth_min_cov: 90.55,
              callable_breadth_min_cov: 100.0,
              mapped_reads: 900,
              qc_status: 'PASS',
              qc_reasons: ['meets_all_pass_thresholds'],
              qc_basis: 'primer_bed',
              coverage_profile: [{ pos: 1, depth: 0 }],
            },
          },
          consensus: {
            FIPV_NC_002306: {
              length: 434,
              gc_content: 44.0,
              n_content: 9.45,
              internal_n_content: 0.0,
              sequence: 'N'.repeat(21) + 'A'.repeat(393) + 'N'.repeat(20),
            },
          },
        },
        {
          sample_id: 'NTC',
          qc_status: 'FAIL',
          qc_reasons: ['callable_breadth_below_warn'],
          primary_reference: 'FIPV_NC_002306',
          qc: {
            raw_reads: 10,
            raw_bases: 3000,
            trimmed_reads: 5,
            trimmed_bases: 1500,
            q30_rate_raw: 0.2,
            q30_rate_trimmed: 0.3,
            duplication_rate: null,
          },
          coverage: {
            FIPV_NC_002306: {
              mean_depth: 0,
              callable_mean_depth: 0,
              coverage_breadth_min_cov: 0,
              callable_breadth_min_cov: 0,
              mapped_reads: 0,
              qc_status: 'FAIL',
              qc_reasons: ['callable_breadth_below_warn'],
              qc_basis: 'primer_bed',
              coverage_profile: [{ pos: 1, depth: 0 }],
            },
          },
          consensus: {
            FIPV_NC_002306: {
              length: 434,
              gc_content: 0,
              n_content: 100,
              internal_n_content: 100,
              sequence: 'N'.repeat(434),
            },
          },
        },
      ],
    };

    expect(() => validateUploadJSON(v2Report)).not.toThrow();
    expect(getSampleStatus(v2Report.samples[0], v2Report).status).toBe('PASS');
    expect(getSampleStatus(v2Report.samples[0], v2Report).isV2).toBe(true);
    expect(getBestReference(v2Report.samples[0], v2Report).name).toBe('FIPV_NC_002306');
    expect(formatMetric(v2Report.samples[0].qc.duplication_rate, { decimals: 1, suffix: '%', isRatioPercent: true })).toBe('N/A');
    expect(formatMetric(v2Report.samples[0].qc.q30_rate_trimmed, { decimals: 1, suffix: '%', isRatioPercent: true })).toBe('56.0%');
  });

  it('falls back cleanly on legacy v1 report without overriding v2 fields', () => {
    const v1Report = {
      run_id: 'legacy_v1',
      min_coverage_threshold: 10,
      references: [{ name: 'SARS2', length: 29903, sequence: 'A'.repeat(100) }],
      samples: [
        {
          sample_id: 'S1',
          qc: {
            raw_reads: 1000,
            raw_bases: 150000,
            trimmed_reads: 950,
            trimmed_bases: 140000,
            q30_rate_raw: 0.9,
            q30_rate_trimmed: 0.92,
            duplication_rate: 0.05,
          },
          coverage: {
            SARS2: {
              mean_depth: 120,
              coverage_breadth_10x: 96.0,
              coverage_profile: [{ pos: 1, depth: 120 }],
            },
          },
          consensus: {
            SARS2: { length: 29903, gc_content: 38.0, n_content: 2.0, sequence: 'A'.repeat(100) },
          },
        },
      ],
    };

    expect(() => validateUploadJSON(v1Report)).not.toThrow();
    const st = getSampleStatus(v1Report.samples[0], v1Report);
    expect(st.status).toBe('PASS');
    expect(st.isV2).toBe(false);
    expect(getMappedReads(v1Report.samples[0], 'SARS2')).toBeNull();
  });

  it('formats overlapping consequences cleanly without p.AA annotation unavailable', () => {
    const variant = {
      gene: 'ORF1ab/ORF1a',
      aa_change: 'S135R',
      consequences: [
        { gene: 'ORF1ab', aa_change: 'AA annotation unavailable', effect: 'unknown', status: 'unsupported_programmed_frameshift' },
        { gene: 'ORF1a', aa_change: 'S135R', effect: 'missense', status: null },
        { gene: 'S', aa_change: 'In-frame del (3 aa)', effect: 'inframe_deletion', status: null },
      ],
    };
    const formatted = formatConsequences(variant);
    expect(formatted[0].label).toBe('AA unavailable (unsupported_programmed_frameshift)');
    expect(formatted[0].kind).toBe('unavailable');
    expect(formatted[1].label).toBe('p.S135R');
    expect(formatted[2].label).toBe('In-frame del (3 aa)');
    expect(formatted[2].kind).toBe('indel');
  });

  it('produces deterministic report IDs and validates bundled mock_run_report.json', () => {
    const id1 = computeStableReportId('runA', 'Sample01');
    const id2 = computeStableReportId('runA', 'Sample01');
    const id3 = computeStableReportId('runB', 'Sample01');
    expect(id1).toBe(id2);
    expect(id1).not.toBe(id3);

    const mockPath = path.resolve(process.cwd(), 'public/mock_run_report.json');
    const mockData = JSON.parse(fs.readFileSync(mockPath, 'utf-8'));
    expect(() => validateUploadJSON(mockData)).not.toThrow();
  });

  it('deduplicates overlapping identical protein consequences and formats depths uniformly', () => {
    const dedup = deduplicateConsequences([
      { gene: 'ORF1ab', product: 'ORF1ab polyprotein', effect: 'missense', status: null, label: 'p.K47R', kind: 'missense' },
      { gene: 'ORF1a', product: 'ORF1a polyprotein', effect: 'missense', status: null, label: 'p.K47R', kind: 'missense' },
      { gene: 'ORF1ab', product: 'ORF1ab polyprotein', effect: 'unknown', status: 'unsupported_programmed_frameshift', label: 'AA unavailable', kind: 'unavailable' },
    ]);
    expect(dedup).toHaveLength(2);
    expect(dedup[0].gene).toBe('ORF1ab/ORF1a');
    expect(dedup[0].label).toBe('p.K47R');
    expect(dedup[1].kind).toBe('unavailable');

    expect(formatDepth(884.18)).toBe('884.2x');
    expect(formatDepth(950.5)).toBe('950.5x');
    expect(formatDepth(null)).toBe('N/A');
  });
});
