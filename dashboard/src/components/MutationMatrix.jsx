import React, { useState } from 'react';
import { Layers, HelpCircle, AlertCircle, ArrowUpDown } from 'lucide-react';

export default function MutationMatrix({ 
  samples, 
  selectedRefName, 
  setSelectedRefName, 
  refs = [],
  minCoverageThreshold = 10
}) {
  const [hoveredCell, setHoveredCell] = useState(null);
  const [sortAsc, setSortAsc] = useState(true);

  if (!samples || samples.length === 0) {
    return null;
  }

  const getBreadth = (sample) => {
    const cov = (selectedRefName && sample.coverage && sample.coverage[selectedRefName]) ? sample.coverage[selectedRefName] : (sample.coverage || {});
    return cov.coverage_breadth_min_cov !== undefined ? cov.coverage_breadth_min_cov : (cov.coverage_breadth_10x !== undefined ? cov.coverage_breadth_10x : 0);
  };

  const getConsensusSeq = (sample, refName) => {
    if (!sample || !sample.consensus) return '';
    const c = (refName && sample.consensus[refName]) ? sample.consensus[refName] : sample.consensus;
    if (!c) return '';
    if (typeof c === 'string') return c;
    if (typeof c.sequence === 'string') return c.sequence;
    return '';
  };

  // Filter samples into passing (>= 50% coverage breadth) and excluded (< 50%) for the active strain
  const passingSamples = samples.filter(sample => getBreadth(sample) >= 50.0);
  const excludedSamples = samples.filter(sample => getBreadth(sample) < 50.0);

  // Sort passing samples alphabetically and numerically
  const sortedPassingSamples = [...passingSamples].sort((a, b) => {
    const comp = a.sample_id.localeCompare(b.sample_id, undefined, { numeric: true, sensitivity: 'base' });
    return sortAsc ? comp : -comp;
  });

  // 1. Collect all unique mutations across passing samples (matching active reference)
  const mutationMap = {};
  passingSamples.forEach(sample => {
    if (sample.variants) {
      sample.variants.forEach(v => {
        if (selectedRefName && v.reference && v.reference !== selectedRefName) return;
        const id = `${v.position}:${v.ref}>${v.alt}`;
        mutationMap[id] = {
          id,
          position: v.position,
          ref: v.ref,
          alt: v.alt,
          gene: v.gene,
          type: v.type,
          aa_change: v.aa_change,
          status: v.status
        };
      });
    }
  });

  // Sort mutations by genomic position
  const sortedMutations = Object.values(mutationMap).sort((a, b) => a.position - b.position);

  // 2. Build sample-to-mutation frequency mapping for quick matrix lookup
  const getCellData = (sampleId, mutId) => {
    const sample = samples.find(s => s.sample_id === sampleId);
    if (!sample || !sample.variants) return null;
    const match = sample.variants.find(v => {
      if (selectedRefName && v.reference && v.reference !== selectedRefName) return false;
      return `${v.position}:${v.ref}>${v.alt}` === mutId;
    });
    return match ? { frequency: match.frequency, depth: match.depth } : null;
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      {/* Header Row with dynamic dropdown */}
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold font-outfit text-white">Mutational Sharing Matrix</h3>
          <p className="text-xs text-gray-400">
            Epidemiological grid comparing variant profiles. Samples with identical columns represent potential transmission clusters.
          </p>
        </div>

        {/* Dynamic Strain Dropdown */}
        {refs && refs.length > 1 && setSelectedRefName && (
          <div className="flex items-center gap-3 shrink-0">
            <span className="text-xs text-gray-400 font-bold uppercase tracking-wider font-outfit">Reference Strain:</span>
            <select
              value={selectedRefName}
              onChange={(e) => setSelectedRefName(e.target.value)}
              className="bg-dark-bg border border-dark-border text-white text-xs font-normal font-outfit rounded-xl px-3 py-1.5 focus:outline-none focus:border-brand-teal cursor-pointer"
            >
              {refs.map(r => (
                <option key={r.name} value={r.name} className="bg-dark-card font-normal">
                  {r.name}
                </option>
              ))}
            </select>
          </div>
        )}
      </div>

      {/* Main Grid Content */}
      {passingSamples.length === 0 ? (
        <div className="p-8 text-center text-gray-400 border border-dark-border rounded-xl bg-dark-bg/20 flex flex-col items-center gap-2">
          <AlertCircle className="w-8 h-8 text-brand-amber animate-pulse" />
          <p className="font-semibold text-white">No samples have sufficient coverage (≥ 50%) for reference: {selectedRefName || 'all'}</p>
          <p className="text-xs text-gray-500">Please select a different reference strain or adjust pipeline filters.</p>
        </div>
      ) : sortedMutations.length === 0 ? (
        <div className="p-8 text-center text-gray-400 border border-dark-border rounded-xl bg-dark-bg/20">
          <p>No variants identified across any of the {passingSamples.length} passing samples for reference: {selectedRefName || 'all'}.</p>
        </div>
      ) : (
        <div className="overflow-x-auto border border-dark-border rounded-xl bg-dark-bg/40 max-w-full">
          <table className="border-collapse w-full min-w-[700px] text-xs">
            <thead>
              <tr className="border-b border-dark-border">
                {/* Top-left cell - Clickable for sorting */}
                <th 
                  className="p-3 text-left font-semibold text-gray-400 border-r border-dark-border w-48 font-outfit shrink-0 cursor-pointer hover:text-white select-none transition-colors"
                  onClick={() => setSortAsc(!sortAsc)}
                >
                  <div className="flex items-center gap-1">
                    Sample ID <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                  </div>
                </th>
                {/* Columns: Mutations */}
                {sortedMutations.map((mut) => (
                  <th 
                    key={mut.id} 
                    className="p-3 text-center font-mono font-medium text-gray-300 border-r border-dark-border min-w-20"
                  >
                    <div className="flex flex-col items-center">
                      <span className="text-[10px] text-brand-indigo font-bold">{mut.gene}</span>
                      <span className="text-white mt-0.5">{mut.ref}{mut.position}{mut.alt}</span>
                      {mut.aa_change && (
                        <span className={`text-[9px] mt-0.5 font-bold ${
                          mut.aa_change.includes('(syn)') 
                            ? 'text-gray-500 font-normal' 
                            : mut.aa_change === 'AA annotation unavailable'
                            ? 'text-gray-500 font-normal italic'
                            : mut.aa_change === 'Indel'
                            ? 'text-brand-amber font-normal'
                            : 'text-brand-teal'
                        }`}>
                          {mut.aa_change.includes('(syn)') 
                            ? mut.aa_change.replace(' (syn)', '') 
                            : mut.aa_change === 'AA annotation unavailable'
                            ? 'unavail'
                            : mut.aa_change}
                        </span>
                      )}
                    </div>
                  </th>
                ))}
              </tr>
            </thead>
            
            <tbody>
              {sortedPassingSamples.map((sample) => (
                <tr 
                  key={sample.sample_id} 
                  className="border-b border-dark-border/40 hover:bg-dark-hover/40 transition-colors"
                >
                  {/* Sample Row Header */}
                  <td className="p-3 font-semibold text-gray-200 border-r border-dark-border font-outfit truncate max-w-48">
                    {sample.sample_id}
                  </td>
                  
                  {/* Cells for mutations */}
                  {sortedMutations.map((mut) => {
                    const cell = getCellData(sample.sample_id, mut.id);
                    const hasMut = cell !== null;
                    const isHighFreq = hasMut && cell.frequency >= 0.90;

                    // Check if consensus base at position is 'N' (low coverage / uncalled)
                    const seq = getConsensusSeq(sample, selectedRefName);
                    const posIdx = mut.position - 1;
                    const baseAtPos = (seq && posIdx >= 0 && posIdx < seq.length && seq[posIdx]) 
                      ? seq[posIdx].toUpperCase() 
                      : null;
                    const isLowCov = !hasMut && (baseAtPos === 'N' || baseAtPos === '-');

                    return (
                      <td 
                        key={mut.id}
                        className="p-1 border-r border-dark-border/40 text-center relative"
                        onMouseEnter={() => setHoveredCell({ sampleId: sample.sample_id, mutId: mut.id, cell, mut, isLowCov })}
                        onMouseLeave={() => setHoveredCell(null)}
                      >
                        <div 
                          className={`mx-auto w-8 h-8 rounded-lg flex items-center justify-center transition-all duration-150 ${
                            hasMut 
                              ? isHighFreq
                                ? 'bg-brand-teal text-dark-bg font-bold shadow-glow-teal scale-90 border border-brand-teal'
                                : 'bg-brand-teal/30 text-brand-teal border border-brand-teal/50 scale-95'
                              : isLowCov
                                ? 'bg-slate-800/80 border border-slate-700/60 text-slate-500 scale-90'
                                : 'bg-white border border-dark-border hover:border-gray-400 hover:bg-dark-hover'
                          }`}
                        >
                          {hasMut ? (
                            <span className="text-[10px] font-mono">
                              {Math.round(cell.frequency * 100)}
                            </span>
                          ) : isLowCov ? (
                            <span className="text-[9px] font-mono text-slate-400 font-bold">
                              N
                            </span>
                          ) : null}
                        </div>
                      </td>
                    );
                  })}
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {/* Grid Legend & Interactive Inspector */}
      {passingSamples.length > 0 && sortedMutations.length > 0 && (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 items-center">
          {/* Cell Details Inspector */}
          <div className="p-4 bg-dark-bg/60 border border-dark-border rounded-xl min-h-20">
            {hoveredCell ? (
              <div className="flex flex-col gap-1">
                <span className="text-[10px] text-gray-400 font-mono">{hoveredCell.sampleId}</span>
                <div className="flex items-center gap-2 mt-0.5">
                  <span className="text-sm font-bold text-white font-mono">
                    g.{hoveredCell.mut.position}{hoveredCell.mut.ref}&gt;{hoveredCell.mut.alt}
                  </span>
                  <span className="text-[10px] bg-brand-teal/20 text-brand-teal font-bold px-1.5 py-0.5 rounded">
                    {hoveredCell.mut.gene}
                  </span>
                  {hoveredCell.mut.aa_change && (
                    <span className={`text-[10px] font-bold px-1.5 py-0.5 rounded ${
                      hoveredCell.mut.aa_change.includes('(syn)') || hoveredCell.mut.aa_change === 'AA annotation unavailable'
                        ? 'bg-gray-800 text-gray-400 border border-dark-border'
                        : hoveredCell.mut.aa_change === 'Indel'
                        ? 'bg-amber-950/40 text-brand-amber border border-brand-amber/30'
                        : 'bg-brand-indigo/20 text-brand-indigo border border-brand-indigo/30'
                    }`}>
                      {hoveredCell.mut.aa_change.includes('(syn)') 
                        ? `${hoveredCell.mut.aa_change.replace(' (syn)', '')} (Synonymous)` 
                        : hoveredCell.mut.aa_change === 'AA annotation unavailable'
                        ? `AA unavailable${hoveredCell.mut.status ? ` (${hoveredCell.mut.status})` : ''}`
                        : hoveredCell.mut.aa_change === 'Indel'
                        ? 'Indel'
                        : `${hoveredCell.mut.aa_change.split('/').map(a => a.startsWith('p.') ? a : `p.${a}`).join('/')} (Non-synonymous)`}
                    </span>
                  )}
                </div>
                {hoveredCell.cell ? (
                  <div className="flex items-center gap-4 text-xs mt-1 text-gray-300">
                    <span>Allele Freq: <strong className="text-brand-teal font-mono">{(hoveredCell.cell.frequency*100).toFixed(1)}%</strong></span>
                    <span>Depth: <strong className="text-white font-mono">{hoveredCell.cell.depth}x</strong></span>
                  </div>
                ) : hoveredCell.isLowCov ? (
                  <div className="text-xs mt-1 text-brand-rose">
                    <span>Low Coverage / Uncalled Region (masked as 'N' in consensus FASTA, &lt;{minCoverageThreshold}x depth)</span>
                  </div>
                ) : (
                  <div className="text-xs mt-1 text-brand-emerald">
                    <span>Wild Type Reference Allele ({hoveredCell.mut.ref})</span>
                  </div>
                )}
              </div>
            ) : (
              <div className="flex items-center gap-2 text-gray-400 text-xs h-full py-2 justify-center">
                <HelpCircle className="w-4 h-4 text-brand-teal" />
                <span>Hover over grid cells to inspect mutation and reference status</span>
              </div>
            )}
          </div>

          {/* Legend */}
          <div className="flex flex-wrap gap-4 justify-start md:justify-end text-xs text-gray-400">
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-brand-teal shadow-glow-teal border border-brand-teal" />
              <span>Clonal Variant (≥90%)</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-brand-teal/30 border border-brand-teal/50" />
              <span>Subclonal Variant (&lt;90%)</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-white border border-dark-border" />
              <span>Wild Type Reference</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-slate-800 border border-slate-700 flex items-center justify-center text-[9px] font-mono text-slate-400 font-bold">N</div>
              <span>Low Coverage / Uncalled (&lt;{minCoverageThreshold}x)</span>
            </div>
          </div>
        </div>
      )}

      {/* Excluded Samples Box */}
      {excludedSamples.length > 0 && (
        <div className="border border-dark-border/60 bg-dark-bg/20 rounded-xl p-4 flex flex-col gap-2">
          <span className="text-xs font-semibold text-brand-rose uppercase tracking-wider font-outfit block">
            Excluded Low-Coverage Samples (&lt;50% Breadth at ≥{minCoverageThreshold}x)
          </span>
          <p className="text-[11px] text-gray-400">
            The following samples had insufficient coverage breadth for this strain and were excluded from comparative mutational sharing analysis to avoid false clustering on uncalled regions (NNNs):
          </p>
          <div className="flex flex-wrap gap-2 mt-1">
            {excludedSamples.map(sample => {
              const breadth = getBreadth(sample);
              return (
                <span key={sample.sample_id} className="text-[11px] bg-dark-card border border-dark-border text-gray-300 font-mono px-2.5 py-1 rounded-lg">
                  {sample.sample_id}: <strong className="text-brand-rose">{breadth.toFixed(1)}% cov</strong>
                </span>
              );
            })}
          </div>
        </div>
      )}
    </div>
  );
}
