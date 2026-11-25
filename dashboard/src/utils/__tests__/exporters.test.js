import { describe, expect, it } from 'vitest';
import {
  escapeCsvField,
  formatCohortMultiFasta,
  formatDistanceMatrixCsv,
  formatSampleFasta,
  formatVariantsCsv,
  treeToNewick,
} from '../exporters.js';

describe('exporters', () => {
  it('escapes CSV fields properly including quotes, commas, and Reference/Consequences columns', () => {
    expect(escapeCsvField('a,b"c')).toBe('"a,b""c"');
    const csv = formatVariantsCsv(
      [
        {
          reference: 'NC_045512.2',
          position: 670,
          ref: 'A',
          alt: 'G',
          gene: 'ORF1ab/ORF1a',
          type: 'Coding',
          depth: 450,
          frequency: 0.99,
          quality: null,
          is_consensus: true,
          aa_change: 'S135R',
          consequences: [{ gene: 'ORF1a', aa_change: 'S135R', effect: 'missense', status: null }],
        },
      ],
      'Sample01'
    );
    const lines = csv.trim().split('\n');
    expect(lines[0]).toContain('Sample_ID,Reference,Position');
    expect(lines[1]).toContain('"Sample01","NC_045512.2",670');
    expect(lines[1]).toContain('"ORF1a:p.S135R"');
  });

  it('formats single and cohort FASTA sequences cleanly with line wrapping', () => {
    const samples = [
      {
        sample_id: 'S1',
        coverage: { ref1: { callable_breadth_min_cov: 95 } },
        consensus: { ref1: { sequence: 'ACGT'.repeat(20) } },
      },
      {
        sample_id: 'S2_low',
        coverage: { ref1: { callable_breadth_min_cov: 40 } },
        consensus: { ref1: { sequence: 'NNNN' } },
      },
    ];
    expect(formatSampleFasta(samples[0], 'ref1').startsWith('>S1_ref1\n')).toBe(true);
    const cohort = formatCohortMultiFasta(samples, 'ref1', { minBreadth: 80 });
    expect(cohort).toContain('>S1_ref1');
    expect(cohort).not.toContain('S2_low');
  });

  it('serializes Newick trees and distance matrices with NA for null entries', () => {
    const tree = {
      id: 'nj_root',
      isLeaf: false,
      children: [
        { id: 'leaf_0', name: 'Sample A', isLeaf: true, branchLength: 0.0123456, children: [] },
        { id: 'leaf_1', name: 'Sample_B', isLeaf: true, branchLength: 0.0123456, children: [] },
      ],
    };
    expect(treeToNewick(tree)).toBe('(Sample_A:0.012346,Sample_B:0.012346);\n');

    const distCsv = formatDistanceMatrixCsv(
      ['A', 'B'],
      [
        [0, null],
        [null, 0],
      ],
      [
        [0, null],
        [null, 0],
      ],
      [
        [100, 0],
        [0, 100],
      ]
    );
    expect(distCsv).toContain('Divergence,"A",0.000000,NA');
    expect(distCsv).toContain('SNP_Diffs,"A",0,NA');
    expect(distCsv).toContain('Compared_Sites,"A",100,0');
  });
});
