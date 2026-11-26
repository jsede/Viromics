import { formatConsequences } from './reportSchema.js';

/**
 * Computes clean power-of-10 logarithmic domain and tick array for depth charts.
 */
export function computeLogScaleConfig(values, threshold = 0) {
  const safeValues = Array.isArray(values)
    ? values.filter((v) => typeof v === 'number' && Number.isFinite(v) && v >= 0)
    : [];
  const safeThreshold = typeof threshold === 'number' && Number.isFinite(threshold) ? threshold : 0;
  const maxVal = Math.max(10, safeThreshold, ...safeValues);
  const maxExp = Math.max(2, Math.ceil(Math.log10(maxVal)));
  const ticks = Array.from({ length: maxExp + 1 }, (_, i) => i);
  return {
    domain: [0, maxExp],
    ticks,
  };
}

/**
 * Formats logarithmic exponents (10^tick) or linear depth ticks into compact human labels.
 */
export function formatDepthAxisTick(tick, isLog = true) {
  if (typeof tick !== 'number' || Number.isNaN(tick)) return '';
  const val = isLog ? Math.round(Math.pow(10, tick)) : tick;
  if (val >= 1_000_000) {
    const m = val / 1_000_000;
    return `${Number.isInteger(m) ? m : m.toFixed(1)}M`;
  }
  if (val >= 1_000) {
    const k = val / 1_000;
    return `${Number.isInteger(k) ? k : k.toFixed(1)}k`;
  }
  return String(Math.round(val));
}

/**
 * Returns Okabe-Ito / Paul Tol derived visual encoding for variant lollipops on GeneMap.
 */
export function getVariantConsequenceStyle(variant) {
  const cons = formatConsequences(variant);
  if (cons.some((c) => c.kind === 'missense' || c.kind === 'nonsense')) {
    return { kind: 'missense', color: '#B43B2E', label: 'Missense / Nonsense' };
  }
  if (cons.some((c) => c.kind === 'frameshift' || c.kind === 'indel')) {
    return { kind: 'indel', color: '#CC8400', label: 'Indel / Frameshift' };
  }
  if (cons.some((c) => c.kind === 'synonymous')) {
    return { kind: 'synonymous', color: '#5C6773', label: 'Synonymous' };
  }
  return { kind: 'noncoding', color: '#6B7280', label: 'Non-coding / Unannotated' };
}
