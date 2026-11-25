import { formatConsequences, getContigConsensus, getEvalBreadth } from './reportSchema.js';

export function escapeCsvField(val) {
  if (val === null || val === undefined) return '""';
  const str = String(val);
  return `"${str.replace(/"/g, '""')}"`;
}

export function formatVariantsCsv(variants, sampleId = '') {
  const headers = [
    'Sample_ID',
    'Reference',
    'Position',
    'Ref',
    'Alt',
    'Gene',
    'Type',
    'Depth',
    'Frequency',
    'Quality',
    'Is_Consensus',
    'AA_Change',
    'Consequences',
    'Status',
  ];
  const rows = [headers.join(',')];

  (variants || []).forEach((v) => {
    const consList = formatConsequences(v);
    const consSummary = consList.map((c) => `${c.gene}:${c.label}`).join('; ');
    rows.push(
      [
        escapeCsvField(sampleId),
        escapeCsvField(v.reference || ''),
        v.position ?? '',
        escapeCsvField(v.ref || ''),
        escapeCsvField(v.alt || ''),
        escapeCsvField(v.gene || ''),
        escapeCsvField(v.type || ''),
        v.depth ?? '',
        v.frequency ?? '',
        v.quality === null || v.quality === undefined ? '' : v.quality,
        v.is_consensus === undefined ? '' : Boolean(v.is_consensus),
        escapeCsvField(v.aa_change || ''),
        escapeCsvField(consSummary),
        escapeCsvField(v.status || ''),
      ].join(',')
    );
  });

  return rows.join('\n') + '\n';
}

export function formatFastaSequence(header, sequence, wrap = 60) {
  const seq = (sequence || '').trim();
  const lines = [`>${header}`];
  for (let i = 0; i < seq.length; i += wrap) {
    lines.push(seq.slice(i, i + wrap));
  }
  return lines.join('\n') + '\n';
}

export function formatSampleFasta(sample, refName) {
  if (!sample) return '';
  const cons = getContigConsensus(sample, refName);
  if (!cons?.sequence) return '';
  const header = refName ? `${sample.sample_id}_${refName}` : sample.sample_id;
  return formatFastaSequence(header, cons.sequence);
}

export function formatCohortMultiFasta(samples, refName, { minBreadth = 0 } = {}) {
  if (!Array.isArray(samples)) return '';
  const chunks = [];
  samples.forEach((sample) => {
    if (minBreadth > 0 && getEvalBreadth(sample, refName) < minBreadth) return;
    const fasta = formatSampleFasta(sample, refName);
    if (fasta) chunks.push(fasta);
  });
  return chunks.join('');
}

function sanitizeNewickLabel(label) {
  return String(label || 'unknown').replace(/[:;,()[\]\s]/g, '_');
}

export function treeToNewick(rootNode) {
  if (!rootNode) return '';
  const serialize = (node, includeBranch) => {
    let sub;
    if (node.isLeaf || !node.children || node.children.length === 0) {
      sub = sanitizeNewickLabel(node.name);
    } else {
      sub = `(${node.children.map((ch) => serialize(ch, true)).join(',')})`;
    }
    if (includeBranch) {
      const len = typeof node.branchLength === 'number' ? Math.max(0, node.branchLength) : 0;
      return `${sub}:${len.toFixed(6)}`;
    }
    return sub;
  };
  return `${serialize(rootNode, false)};\n`;
}

export function formatDistanceMatrixCsv(sampleIds, distMatrix, snpMatrix = null, sitesMatrix = null) {
  const n = sampleIds.length;
  const rows = [];
  rows.push(['Metric', 'Sample_ID', ...sampleIds.map(escapeCsvField)].join(','));

  for (let i = 0; i < n; i++) {
    const vals = [];
    for (let j = 0; j < n; j++) {
      const v = distMatrix?.[i]?.[j];
      vals.push(v === null || v === undefined ? 'NA' : Number(v).toFixed(6));
    }
    rows.push(['Divergence', escapeCsvField(sampleIds[i]), ...vals].join(','));
  }

  if (Array.isArray(snpMatrix) && snpMatrix.length === n) {
    for (let i = 0; i < n; i++) {
      const vals = [];
      for (let j = 0; j < n; j++) {
        const v = snpMatrix?.[i]?.[j];
        vals.push(v === null || v === undefined ? 'NA' : String(v));
      }
      rows.push(['SNP_Diffs', escapeCsvField(sampleIds[i]), ...vals].join(','));
    }
  }

  if (Array.isArray(sitesMatrix) && sitesMatrix.length === n) {
    for (let i = 0; i < n; i++) {
      const vals = [];
      for (let j = 0; j < n; j++) {
        const v = sitesMatrix?.[i]?.[j];
        vals.push(v === null || v === undefined ? 'NA' : String(v));
      }
      rows.push(['Compared_Sites', escapeCsvField(sampleIds[i]), ...vals].join(','));
    }
  }

  return rows.join('\n') + '\n';
}

export function serializeSvgElement(svgElement, { background = '#FFFFFF' } = {}) {
  if (!svgElement || typeof document === 'undefined') return '';
  const clone = svgElement.cloneNode(true);
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.style.fontFamily = 'Arial, Helvetica, sans-serif';

  const vb = clone.getAttribute('viewBox');
  if (vb && background) {
    const parts = vb.trim().split(/\s+/).map(Number);
    if (parts.length === 4 && parts.every((x) => !Number.isNaN(x))) {
      const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      bg.setAttribute('x', String(parts[0]));
      bg.setAttribute('y', String(parts[1]));
      bg.setAttribute('width', String(parts[2]));
      bg.setAttribute('height', String(parts[3]));
      bg.setAttribute('fill', background);
      clone.insertBefore(bg, clone.firstChild);
    }
  }
  return `<?xml version="1.0" encoding="UTF-8"?>\n${clone.outerHTML}\n`;
}

export function triggerBrowserDownload(content, filename, mimeType = 'text/plain;charset=utf-8;') {
  if (typeof document === 'undefined' || typeof URL === 'undefined') return;
  const blob = new Blob([content], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.setAttribute('href', url);
  link.setAttribute('download', filename);
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
