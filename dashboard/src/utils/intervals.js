import { getContigCoverage, getContigConsensus } from './reportSchema.js';

/**
 * Binary search on sorted 1-based inclusive [start, end] intervals.
 */
export function isPositionInIntervals(pos, intervals) {
  if (!Array.isArray(intervals) || intervals.length === 0 || typeof pos !== 'number') {
    return false;
  }
  let lo = 0;
  let hi = intervals.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >> 1;
    const [start, end] = intervals[mid];
    if (pos < start) {
      hi = mid - 1;
    } else if (pos > end) {
      lo = mid + 1;
    } else {
      return true;
    }
  }
  return false;
}

/**
 * Indel-safe low-coverage check at 1-based reference position.
 * Uses v2 reference-coordinate `low_coverage_intervals` when present;
 * falls back to consensus string indexing ONLY on legacy v1 reports.
 */
export function isVariantPositionLowCoverage(sample, refName, position) {
  const cov = getContigCoverage(sample, refName);
  if (Array.isArray(cov.low_coverage_intervals)) {
    return {
      isLowCov: isPositionInIntervals(position, cov.low_coverage_intervals),
      usedFallback: false,
    };
  }

  const cons = getContigConsensus(sample, refName);
  const seq = typeof cons.sequence === 'string' ? cons.sequence : '';
  const idx = position - 1;
  const ch = seq && idx >= 0 && idx < seq.length ? seq[idx].toUpperCase() : null;
  return {
    isLowCov: ch === 'N' || ch === '-',
    usedFallback: true,
  };
}

/**
 * Compute internal low-coverage dropout intervals within callable target regions.
 * No hardcoded terminal buffer: relies directly on pipeline callable_target_intervals.
 */
export function getInternalDropoutIntervals(
  lowCoverageIntervals,
  callableTargetIntervals,
  coverageProfile,
  minCoverageThreshold = 10
) {
  if (Array.isArray(lowCoverageIntervals)) {
    if (Array.isArray(callableTargetIntervals) && callableTargetIntervals.length > 0) {
      const intersections = [];
      lowCoverageIntervals.forEach(([lStart, lEnd]) => {
        callableTargetIntervals.forEach(([cStart, cEnd]) => {
          const s = Math.max(lStart, cStart);
          const e = Math.min(lEnd, cEnd);
          if (s <= e) {
            intersections.push({ start: s, end: e });
          }
        });
      });
      return intersections;
    }
    return lowCoverageIntervals.map(([start, end]) => ({ start, end }));
  }

  // Legacy v1 fallback from binned profile
  if (!Array.isArray(coverageProfile) || coverageProfile.length === 0) {
    return [];
  }
  const dropouts = [];
  let current = null;
  coverageProfile.forEach((pt) => {
    if (pt.depth < minCoverageThreshold) {
      if (!current) {
        current = { start: pt.pos, end: pt.pos };
      } else {
        current.end = pt.pos;
      }
    } else if (current) {
      dropouts.push(current);
      current = null;
    }
  });
  if (current) dropouts.push(current);
  return dropouts;
}
