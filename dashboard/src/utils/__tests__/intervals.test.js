import { describe, expect, it } from 'vitest';
import {
  getInternalDropoutIntervals,
  isPositionInIntervals,
  isVariantPositionLowCoverage,
} from '../intervals.js';

describe('intervals (indel-safe coverage lookup)', () => {
  it('performs binary search over sorted 1-based inclusive intervals', () => {
    const intervals = [
      [1, 54],
      [11288, 11296],
      [29837, 29870],
    ];
    expect(isPositionInIntervals(1, intervals)).toBe(true);
    expect(isPositionInIntervals(54, intervals)).toBe(true);
    expect(isPositionInIntervals(55, intervals)).toBe(false);
    expect(isPositionInIntervals(11290, intervals)).toBe(true);
    expect(isPositionInIntervals(29836, intervals)).toBe(false);
    expect(isPositionInIntervals(29837, intervals)).toBe(true);
  });

  it('prevents downstream SNP shift after deletions by using v2 low_coverage_intervals', () => {
    // Suppose sample has a 55 bp upstream deletion, so consensus length is 200 bp instead of 255 bp,
    // and consensus[240 - 1] is out-of-bounds or shifted to a neighboring 'N'.
    // v2 low_coverage_intervals accurately reports pos 240 as well-covered!
    const sampleV2 = {
      coverage: {
        'NC_045512.2': {
          low_coverage_intervals: [
            [1, 54],
            [250, 255],
          ],
        },
      },
      consensus: {
        'NC_045512.2': {
          sequence: 'A'.repeat(200),
        },
      },
    };

    const resCovered = isVariantPositionLowCoverage(sampleV2, 'NC_045512.2', 240);
    expect(resCovered).toEqual({ isLowCov: false, usedFallback: false });

    const resMasked = isVariantPositionLowCoverage(sampleV2, 'NC_045512.2', 252);
    expect(resMasked).toEqual({ isLowCov: true, usedFallback: false });

    // Legacy v1 without low_coverage_intervals falls back cleanly
    const sampleV1 = {
      coverage: { ref1: { mean_depth: 100 } },
      consensus: { ref1: { sequence: 'AANAA' } },
    };
    expect(isVariantPositionLowCoverage(sampleV1, 'ref1', 3)).toEqual({
      isLowCov: true,
      usedFallback: true,
    });
    expect(isVariantPositionLowCoverage(sampleV1, 'ref1', 2)).toEqual({
      isLowCov: false,
      usedFallback: true,
    });
  });

  it('intersects low_coverage_intervals strictly with callable_target_intervals (ignoring terminal primer ends)', () => {
    // Pan-CoV 434-bp contig: terminal primers [1..21] and [415..434] excluded from callable target [22..414]
    const lowCov = [
      [1, 21],
      [150, 175],
      [415, 434],
    ];
    const callable = [[22, 414]];
    expect(getInternalDropoutIntervals(lowCov, callable, [], 50)).toEqual([
      { start: 150, end: 175 },
    ]);

    // When only terminal primer ends are low coverage, internal dropouts = []
    expect(getInternalDropoutIntervals([[1, 21], [415, 434]], callable, [], 50)).toEqual([]);
  });
});
