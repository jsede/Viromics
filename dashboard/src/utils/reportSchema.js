/**
 * Schema v2.0.0 report accessors, v1 graceful fallbacks, validation, and formatting helpers.
 */

export const DEFAULT_QC_THRESHOLDS = {
  pass_depth: 100,
  pass_breadth: 95.0,
  pass_max_n: 5.0,
  warn_depth: 25,
  warn_breadth: 80.0,
  warn_max_n: 10.0,
};

const VALID_QC_TIERS = new Set(['PASS', 'WARN', 'FAIL']);

function isNumOrNull(val) {
  return val === null || (typeof val === 'number' && !Number.isNaN(val));
}

export function validateUploadJSON(data) {
  if (!data || typeof data !== 'object' || Array.isArray(data)) {
    throw new Error('Report must be a JSON object.');
  }
  if (!data.run_id || typeof data.run_id !== 'string') {
    throw new Error("Missing or invalid 'run_id' in report.");
  }
  if (data.min_coverage_threshold === undefined || typeof data.min_coverage_threshold !== 'number') {
    throw new Error("Missing or invalid 'min_coverage_threshold' in report.");
  }

  const refs = getReferencesList(data);
  if (refs.length === 0) {
    throw new Error("Missing or invalid 'references' list in report.");
  }
  refs.forEach((ref, idx) => {
    if (!ref.name || typeof ref.name !== 'string') {
      throw new Error(`Reference at index ${idx} is missing 'name'.`);
    }
    if (typeof ref.length !== 'number' || ref.length < 0) {
      throw new Error(`Reference '${ref.name}' is missing valid 'length'.`);
    }
  });

  if (!data.samples || !Array.isArray(data.samples)) {
    throw new Error("Missing or invalid 'samples' list in report.");
  }

  data.samples.forEach((s, idx) => {
    const sId = s.sample_id || `index_${idx}`;
    if (!s.sample_id || typeof s.sample_id !== 'string') {
      throw new Error(`Sample at index ${idx} is missing 'sample_id'.`);
    }
    if (!s.qc || typeof s.qc !== 'object') {
      throw new Error(`Sample '${sId}' is missing 'qc' dictionary.`);
    }
    if (!s.coverage || typeof s.coverage !== 'object') {
      throw new Error(`Sample '${sId}' is missing 'coverage' dictionary.`);
    }
    if (!s.consensus || typeof s.consensus !== 'object') {
      throw new Error(`Sample '${sId}' is missing 'consensus' dictionary.`);
    }

    ['raw_reads', 'raw_bases', 'trimmed_reads', 'trimmed_bases'].forEach((key) => {
      if (typeof s.qc[key] !== 'number' || Number.isNaN(s.qc[key])) {
        throw new Error(`Sample '${sId}' QC is missing or has invalid number for '${key}'.`);
      }
    });

    ['q30_rate_raw', 'q30_rate_trimmed', 'duplication_rate'].forEach((key) => {
      if (s.qc[key] === undefined || !isNumOrNull(s.qc[key])) {
        throw new Error(`Sample '${sId}' QC has invalid value for nullable metric '${key}'.`);
      }
    });

    const covEntries =
      s.coverage.mean_depth !== undefined ? { __flat__: s.coverage } : s.coverage;
    Object.entries(covEntries).forEach(([refName, cov]) => {
      if (!cov || typeof cov !== 'object') {
        throw new Error(`Sample '${sId}' is missing coverage data for reference '${refName}'.`);
      }
      if (typeof cov.mean_depth !== 'number') {
        throw new Error(`Sample '${sId}' coverage for '${refName}' is missing 'mean_depth'.`);
      }
      if (!Array.isArray(cov.coverage_profile)) {
        throw new Error(`Sample '${sId}' coverage_profile for '${refName}' must be a list.`);
      }
      if (
        cov.coverage_breadth_min_cov === undefined &&
        cov.coverage_breadth_10x === undefined &&
        cov.callable_breadth_min_cov === undefined
      ) {
        throw new Error(`Sample '${sId}' coverage for '${refName}' is missing breadth metrics.`);
      }
    });

    const consEntries =
      s.consensus.sequence !== undefined ? { __flat__: s.consensus } : s.consensus;
    Object.entries(consEntries).forEach(([refName, cons]) => {
      if (!cons || typeof cons !== 'object') {
        throw new Error(`Sample '${sId}' is missing consensus data for reference '${refName}'.`);
      }
      ['length', 'gc_content', 'n_content', 'sequence'].forEach((key) => {
        if (cons[key] === undefined) {
          throw new Error(`Sample '${sId}' consensus for '${refName}' is missing '${key}'.`);
        }
      });
    });
  });

  if (data.distance_matrices && typeof data.distance_matrices === 'object') {
    Object.entries(data.distance_matrices).forEach(([refName, rows]) => {
      if (!Array.isArray(rows)) {
        throw new Error(`Invalid distance_matrices for '${refName}'.`);
      }
      rows.forEach((row) => {
        if (!Array.isArray(row) || row.some((v) => !isNumOrNull(v))) {
          throw new Error(`Invalid distance entry in distance_matrices['${refName}'].`);
        }
      });
    });
  }
}

export function getReferencesList(report) {
  if (!report) return [];
  if (Array.isArray(report.references)) return report.references;
  if (report.reference) return [report.reference];
  return [];
}

export function getRunConfig(report) {
  return report?.run_config && typeof report.run_config === 'object' ? report.run_config : {};
}

export function getQcThresholds(report) {
  const raw = report?.run_config?.qc_thresholds;
  if (!raw || typeof raw !== 'object') return { ...DEFAULT_QC_THRESHOLDS };
  return {
    pass_depth: raw.pass_depth ?? DEFAULT_QC_THRESHOLDS.pass_depth,
    pass_breadth: raw.pass_breadth ?? DEFAULT_QC_THRESHOLDS.pass_breadth,
    pass_max_n: raw.pass_max_n ?? DEFAULT_QC_THRESHOLDS.pass_max_n,
    warn_depth: raw.warn_depth ?? DEFAULT_QC_THRESHOLDS.warn_depth,
    warn_breadth: raw.warn_breadth ?? DEFAULT_QC_THRESHOLDS.warn_breadth,
    warn_max_n: raw.warn_max_n ?? DEFAULT_QC_THRESHOLDS.warn_max_n,
  };
}

export function getSoftwareVersions(report) {
  return report?.software_versions && typeof report.software_versions === 'object'
    ? report.software_versions
    : {};
}

export function getContigCoverage(sample, refName) {
  if (!sample?.coverage) return {};
  if (refName && sample.coverage[refName]) return sample.coverage[refName];
  if (sample.coverage.mean_depth !== undefined) return sample.coverage;
  return {};
}

export function getContigConsensus(sample, refName) {
  if (!sample?.consensus) return {};
  if (refName && sample.consensus[refName]) return sample.consensus[refName];
  if (sample.consensus.sequence !== undefined) return sample.consensus;
  if (typeof sample.consensus === 'string') {
    return { sequence: sample.consensus, length: sample.consensus.length };
  }
  return {};
}

export function getEvalBreadth(sample, refName) {
  const cov = getContigCoverage(sample, refName);
  if (typeof cov.callable_breadth_min_cov === 'number') return cov.callable_breadth_min_cov;
  if (typeof cov.coverage_breadth_min_cov === 'number') return cov.coverage_breadth_min_cov;
  if (typeof cov.coverage_breadth_10x === 'number') return cov.coverage_breadth_10x;
  return 0;
}

export function getWholeBreadth(sample, refName) {
  const cov = getContigCoverage(sample, refName);
  if (typeof cov.coverage_breadth_min_cov === 'number') return cov.coverage_breadth_min_cov;
  if (typeof cov.coverage_breadth_10x === 'number') return cov.coverage_breadth_10x;
  return 0;
}

export function getEvalDepth(sample, refName) {
  const cov = getContigCoverage(sample, refName);
  if (typeof cov.callable_mean_depth === 'number') return cov.callable_mean_depth;
  if (typeof cov.mean_depth === 'number') return cov.mean_depth;
  return 0;
}

export function getEvalNRate(sample, refName) {
  const cons = getContigConsensus(sample, refName);
  if (typeof cons.internal_n_content === 'number') return cons.internal_n_content;
  if (typeof cons.n_content === 'number') return cons.n_content;
  return 0;
}

export function getMappedReads(sample, refName) {
  const cov = getContigCoverage(sample, refName);
  return typeof cov.mapped_reads === 'number' ? cov.mapped_reads : null;
}

export function getBestReference(sample, report) {
  const refs = getReferencesList(report);
  if (refs.length === 0) return null;
  if (sample?.primary_reference) {
    const exact = refs.find((r) => r.name === sample.primary_reference);
    if (exact) return exact;
  }
  if (refs.length === 1) return refs[0];

  let bestRef = refs[0];
  let bestBreadth = -1;
  let bestDepth = -1;
  refs.forEach((r) => {
    const br = getEvalBreadth(sample, r.name);
    const dp = getEvalDepth(sample, r.name);
    if (br > bestBreadth || (br === bestBreadth && dp > bestDepth)) {
      bestBreadth = br;
      bestDepth = dp;
      bestRef = r;
    }
  });
  return bestRef;
}

export function getSampleStatusForRef(sample, refName, report) {
  const cov = getContigCoverage(sample, refName);
  if (cov && VALID_QC_TIERS.has(cov.qc_status)) {
    return {
      status: cov.qc_status,
      reasons: Array.isArray(cov.qc_reasons) ? cov.qc_reasons : [],
      basis: cov.qc_basis || 'whole_contig',
      isV2: true,
    };
  }

  const t = getQcThresholds(report);
  const depth = getEvalDepth(sample, refName);
  const breadth = getEvalBreadth(sample, refName);
  const nRate = getEvalNRate(sample, refName);

  if (depth >= t.pass_depth && breadth >= t.pass_breadth && nRate <= t.pass_max_n) {
    return { status: 'PASS', reasons: ['meets_all_pass_thresholds'], basis: 'whole_contig', isV2: false };
  }
  const failReasons = [];
  if (depth < t.warn_depth) failReasons.push(`mean_depth_below_warn (${depth.toFixed(1)} < ${t.warn_depth})`);
  if (breadth < t.warn_breadth) failReasons.push(`breadth_below_warn (${breadth.toFixed(1)}% < ${t.warn_breadth}%)`);
  if (nRate > t.warn_max_n) failReasons.push(`n_rate_above_warn (${nRate.toFixed(1)}% > ${t.warn_max_n}%)`);
  if (failReasons.length > 0) {
    return { status: 'FAIL', reasons: failReasons, basis: 'whole_contig', isV2: false };
  }
  const warnReasons = [];
  if (depth < t.pass_depth) warnReasons.push(`mean_depth_below_pass (${depth.toFixed(1)} < ${t.pass_depth})`);
  if (breadth < t.pass_breadth) warnReasons.push(`breadth_below_pass (${breadth.toFixed(1)}% < ${t.pass_breadth}%)`);
  if (nRate > t.pass_max_n) warnReasons.push(`n_rate_above_pass (${nRate.toFixed(1)}% > ${t.pass_max_n}%)`);
  return { status: 'WARN', reasons: warnReasons, basis: 'whole_contig', isV2: false };
}

export function getSampleStatus(sample, report) {
  const bestRef = getBestReference(sample, report);
  if (sample && VALID_QC_TIERS.has(sample.qc_status)) {
    const cov = bestRef ? getContigCoverage(sample, bestRef.name) : {};
    return {
      status: sample.qc_status,
      reasons: Array.isArray(sample.qc_reasons) ? sample.qc_reasons : [],
      primaryReference: sample.primary_reference || bestRef?.name || null,
      basis: cov.qc_basis || 'whole_contig',
      isV2: true,
    };
  }
  if (!bestRef) {
    return {
      status: 'FAIL',
      reasons: ['no_reference_coverage'],
      primaryReference: null,
      basis: 'whole_contig',
      isV2: false,
    };
  }
  const res = getSampleStatusForRef(sample, bestRef.name, report);
  return {
    ...res,
    primaryReference: bestRef.name,
  };
}

export function getCoDetections(sample, report, minBreadth = 10.0) {
  const refs = getReferencesList(report);
  if (refs.length <= 1) return [];
  return refs.filter((r) => getEvalBreadth(sample, r.name) >= minBreadth).map((r) => r.name);
}

export function formatMetric(val, { decimals = 2, suffix = '', isRatioPercent = false } = {}) {
  if (val === null || val === undefined || typeof val === 'boolean') return 'N/A';
  const num = Number(val);
  if (Number.isNaN(num)) return 'N/A';
  const scaled = isRatioPercent && num <= 1.0 ? num * 100 : num;
  return `${scaled.toFixed(decimals)}${suffix}`;
}

export function computeStableReportId(runId, sampleId) {
  const input = `${runId || 'run'}::${sampleId || 'sample'}`;
  let hash = 0x811c9dc5;
  for (let i = 0; i < input.length; i++) {
    hash ^= input.charCodeAt(i);
    hash = Math.imul(hash, 0x01000193);
  }
  const hex = (hash >>> 0).toString(16).toUpperCase().padStart(8, '0').slice(0, 6);
  return `VMR-${sampleId || 'UNKNOWN'}-${hex}`;
}

function classifySingleAaLabel(aaChange, effect, status) {
  if (!aaChange) {
    if (status) {
      return { label: `AA unavailable (${status})`, kind: 'unavailable' };
    }
    return { label: '-', kind: 'none' };
  }
  if (aaChange === 'AA annotation unavailable') {
    return {
      label: `AA unavailable${status ? ` (${status})` : ''}`,
      kind: 'unavailable',
    };
  }
  if (
    aaChange === 'Indel' ||
    aaChange.startsWith('In-frame ') ||
    aaChange.startsWith('Frameshift ') ||
    aaChange.startsWith('MNP ')
  ) {
    const isFrameshift = aaChange.startsWith('Frameshift ') || (effect && effect.includes('frameshift'));
    return { label: aaChange, kind: isFrameshift ? 'frameshift' : 'indel' };
  }
  if (aaChange.includes('(syn)') || effect === 'synonymous') {
    return {
      label: `${aaChange.replace(' (syn)', '')} (syn)`,
      kind: 'synonymous',
    };
  }
  const formatted = aaChange
    .split('/')
    .map((part) => (part.startsWith('p.') ? part : `p.${part}`))
    .join('/');
  return { label: formatted, kind: effect === 'nonsense' ? 'nonsense' : 'missense' };
}

export function formatConsequences(variant) {
  if (!variant) return [];
  if (Array.isArray(variant.consequences) && variant.consequences.length > 0) {
    return variant.consequences.map((c) => {
      const { label, kind } = classifySingleAaLabel(c.aa_change, c.effect, c.status);
      return {
        gene: c.gene || variant.gene || 'Unknown',
        proteinId: c.protein_id || null,
        product: c.product || null,
        aaChange: c.aa_change || '',
        effect: c.effect || 'unknown',
        status: c.status || null,
        label,
        kind,
      };
    });
  }

  if (!variant.gene || variant.gene === 'Intergenic') {
    return [];
  }

  const { label, kind } = classifySingleAaLabel(variant.aa_change, null, variant.status);
  return [
    {
      gene: variant.gene,
      proteinId: null,
      product: null,
      aaChange: variant.aa_change || '',
      effect: 'unknown',
      status: variant.status || null,
      label,
      kind,
    },
  ];
}
