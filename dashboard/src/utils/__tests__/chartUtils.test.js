import { describe, expect, it } from 'vitest';
import { computeLogScaleConfig, formatDepthAxisTick, getVariantConsequenceStyle } from '../chartUtils.js';

describe('chartUtils', () => {
  it('computes clean decade domains and ticks across depth scales', () => {
    expect(computeLogScaleConfig([0, 12, 850], 10)).toEqual({
      domain: [0, 3],
      ticks: [0, 1, 2, 3],
    });
    expect(computeLogScaleConfig([100, 12400], 10)).toEqual({
      domain: [0, 5],
      ticks: [0, 1, 2, 3, 4, 5],
    });
    expect(computeLogScaleConfig([], 10)).toEqual({
      domain: [0, 2],
      ticks: [0, 1, 2],
    });
  });

  it('formats logarithmic and linear depth ticks into compact human labels', () => {
    expect([0, 1, 2, 3, 4].map((t) => formatDepthAxisTick(t, true))).toEqual([
      '1',
      '10',
      '100',
      '1k',
      '10k',
    ]);
    expect(formatDepthAxisTick(6, true)).toBe('1M');
    expect(formatDepthAxisTick(2500, false)).toBe('2.5k');
  });

  it('classifies variant lollipop consequence colors consistently', () => {
    expect(
      getVariantConsequenceStyle({
        gene: 'S',
        consequences: [{ gene: 'S', aa_change: 'D614G', effect: 'missense', status: null }],
      }).color
    ).toBe('#B43B2E');
    expect(
      getVariantConsequenceStyle({
        gene: 'S',
        consequences: [{ gene: 'S', aa_change: 'Frameshift ins', effect: 'frameshift', status: null }],
      }).color
    ).toBe('#CC8400');
    expect(
      getVariantConsequenceStyle({
        gene: 'N',
        consequences: [{ gene: 'N', aa_change: 'P13P (syn)', effect: 'synonymous', status: null }],
      }).color
    ).toBe('#5C6773');
    expect(getVariantConsequenceStyle({ gene: 'Intergenic', aa_change: '' }).color).toBe('#6B7280');
  });
});
