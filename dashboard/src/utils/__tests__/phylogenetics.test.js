import { describe, expect, it } from 'vitest';
import {
  buildNeighborJoiningTree,
  buildUPGMATree,
  filterTreeSamples,
} from '../phylogenetics.js';

function getLeafDepths(node, out = {}) {
  if (!node) return out;
  if (node.isLeaf) {
    out[node.name] = node.depthFromRoot;
  } else {
    node.children.forEach((c) => getLeafDepths(c, out));
  }
  return out;
}

function getAllBranches(node, out = []) {
  if (!node) return out;
  out.push(node.branchLength);
  (node.children || []).forEach((c) => getAllBranches(c, out));
  return out;
}

describe('phylogenetics', () => {
  it('builds exact midpoint-rooted Neighbor-Joining tree and clamps negative branches', () => {
    // True additive 4-taxon unrooted tree:
    // (A:0.01, B:0.03) --0.02-- (C:0.04, D:0.05)
    // Longest leaf-to-leaf path: B <-> D = 0.03 + 0.02 + 0.05 = 0.10 (midpoint = 0.05)
    const ids = ['SampleA', 'SampleB', 'SampleC', 'SampleD'];
    const D = [
      [0.0, 0.04, 0.07, 0.08],
      [0.04, 0.0, 0.09, 0.10],
      [0.07, 0.09, 0.0, 0.09],
      [0.08, 0.10, 0.09, 0.0],
    ];

    const nj = buildNeighborJoiningTree(ids, D);
    expect(nj.id).toBe('nj_root');
    const depths = getLeafDepths(nj);
    expect(depths.SampleB).toBeCloseTo(0.05, 6);
    expect(depths.SampleD).toBeCloseTo(0.05, 6);
    expect(depths.SampleA).toBeCloseTo(0.03, 6);
    expect(depths.SampleC).toBeCloseTo(0.04, 6);
    getAllBranches(nj).forEach((len) => expect(len).toBeGreaterThanOrEqual(0));

    // Negative branch clamping test on slightly non-additive distances
    const nonAdditive = [
      [0.0, 0.001, 0.05, 0.06],
      [0.001, 0.0, 0.08, 0.09],
      [0.05, 0.08, 0.0, 0.01],
      [0.06, 0.09, 0.01, 0.0],
    ];
    const njClamped = buildNeighborJoiningTree(ids, nonAdditive);
    getAllBranches(njClamped).forEach((len) => expect(len).toBeGreaterThanOrEqual(0));
  });

  it('builds ultrametric UPGMA tree cleanly', () => {
    const ids = ['S1', 'S2', 'S3'];
    const D = [
      [0.0, 0.02, 0.06],
      [0.02, 0.0, 0.06],
      [0.06, 0.06, 0.0],
    ];
    const upgma = buildUPGMATree(ids, D);
    const depths = getLeafDepths(upgma);
    expect(depths.S1).toBeCloseTo(0.03, 6);
    expect(depths.S2).toBeCloseTo(0.03, 6);
    expect(depths.S3).toBeCloseTo(0.03, 6);
  });

  it('filters low-breadth samples, null pairwise distances, and insufficient compared_sites', () => {
    const samples = [
      { sample_id: 'A', coverage: { ref1: { callable_breadth_min_cov: 98 } } },
      { sample_id: 'B', coverage: { ref1: { callable_breadth_min_cov: 96 } } },
      { sample_id: 'NoOverlap', coverage: { ref1: { callable_breadth_min_cov: 85 } } },
      { sample_id: 'LowBreadth', coverage: { ref1: { callable_breadth_min_cov: 45 } } },
    ];
    const report = {
      distance_matrix_order: ['A', 'B', 'NoOverlap', 'LowBreadth'],
      distance_matrices: {
        ref1: [
          [0.0, 0.001, null, null],
          [0.001, 0.0, null, null],
          [null, null, 0.0, null],
          [null, null, null, 0.0],
        ],
      },
      snp_distance_matrices: {
        ref1: [
          [0, 10, null, null],
          [10, 0, null, null],
          [null, null, 0, null],
          [null, null, null, 0],
        ],
      },
      compared_sites_matrices: {
        ref1: [
          [10000, 9950, 0, 0],
          [9950, 10000, 0, 0],
          [0, 0, 8500, 0],
          [0, 0, 0, 4500],
        ],
      },
    };

    const res = filterTreeSamples(samples, 'ref1', report, { minBreadth: 80, minComparedSites: 100 });
    expect(res.includedSamples.map((s) => s.sample_id)).toEqual(['A', 'B']);
    expect(res.excludedLowBreadth.map((e) => e.sample.sample_id)).toEqual(['LowBreadth']);
    expect(res.excludedOverlap.map((e) => e.sample.sample_id)).toEqual(['NoOverlap']);
    expect(res.distMatrix).toEqual([
      [0, 0.001],
      [0.001, 0],
    ]);
    expect(res.snpMatrix[0][1]).toBe(10);
    expect(res.sitesMatrix[0][1]).toBe(9950);
  });
});
