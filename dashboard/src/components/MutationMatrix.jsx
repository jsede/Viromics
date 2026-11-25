import React, { useMemo, useState } from 'react';
import { HelpCircle, AlertCircle, ArrowUpDown, SlidersHorizontal } from 'lucide-react';
import { formatConsequences, getEvalBreadth } from '../utils/reportSchema';
import { isVariantPositionLowCoverage } from '../utils/intervals';

export default function MutationMatrix({
  samples,
  selectedRefName,
  setSelectedRefName,
  refs = [],
  minCoverageThreshold = 10,
  runConfig = {},
}) {
  const defaultMinFreq = typeof runConfig.min_freq === 'number' ? runConfig.min_freq : 0.2;
  const defaultConsFreq = typeof runConfig.min_consensus_freq === 'number' ? runConfig.min_consensus_freq : 0.5;

  const [hoveredCell, setHoveredCell] = useState(null);
  const [sortAsc, setSortAsc] = useState(true);
  const [minBreadthCutoff, setMinBreadthCutoff] = useState(50);
  const [clonalAfPct, setClonalAfPct] = useState(90);
  const [variantScope, setVariantScope] = useState('ALL'); // 'ALL' | 'CONSENSUS'

  const clonalAf = Math.max(0.01, Math.min(1.0, clonalAfPct / 100));

  const { passingSamples, excludedSamples, sortedPassingSamples } = useMemo(() => {
    const safe = Array.isArray(samples) ? samples : [];
    const passing = [];
    const excluded = [];
    safe.forEach((s) => {
      const br = getEvalBreadth(s, selectedRefName);
      if (br >= minBreadthCutoff) {
        passing.push(s);
      } else {
        excluded.push({ sample: s, breadth: br });
      }
    });
    const sorted = [...passing].sort((a, b) => {
      const comp = a.sample_id.localeCompare(b.sample_id, undefined, { numeric: true, sensitivity: 'base' });
      return sortAsc ? comp : -comp;
    });
    return { passingSamples: passing, excludedSamples: excluded, sortedPassingSamples: sorted };
  }, [samples, selectedRefName, minBreadthCutoff, sortAsc]);

  const sortedMutations = useMemo(() => {
    const map = new Map();
    passingSamples.forEach((sample) => {
      (sample.variants || []).forEach((v) => {
        if (selectedRefName && v.reference && v.reference !== selectedRefName) return;
        const isCons = v.is_consensus !== undefined ? Boolean(v.is_consensus) : v.frequency >= defaultConsFreq;
        if (variantScope === 'CONSENSUS' && !isCons) return;
        const id = `${v.position}:${v.ref}>${v.alt}`;
        if (!map.has(id)) {
          map.set(id, {
            id,
            position: v.position,
            ref: v.ref,
            alt: v.alt,
            gene: v.gene,
            type: v.type,
            aa_change: v.aa_change,
            status: v.status,
            consequences: formatConsequences(v),
          });
        }
      });
    });
    return [...map.values()].sort((a, b) => a.position - b.position);
  }, [passingSamples, selectedRefName, variantScope, defaultConsFreq]);

  const variantLookup = useMemo(() => {
    const lookup = new Map();
    (samples || []).forEach((s) => {
      const perSample = new Map();
      (s.variants || []).forEach((v) => {
        if (selectedRefName && v.reference && v.reference !== selectedRefName) return;
        perSample.set(`${v.position}:${v.ref}>${v.alt}`, {
          frequency: v.frequency,
          depth: v.depth,
          is_consensus: v.is_consensus !== undefined ? Boolean(v.is_consensus) : v.frequency >= defaultConsFreq,
        });
      });
      lookup.set(s.sample_id, perSample);
    });
    return lookup;
  }, [samples, selectedRefName, defaultConsFreq]);

  if (!samples || samples.length === 0) {
    return null;
  }

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold text-slate-900">Mutational Sharing Matrix</h3>
          <p className="text-xs text-slate-500">
            Reference-coordinate mutation matrix using indel-safe low-coverage interval lookup (reported AF &ge; {(defaultMinFreq * 100).toFixed(0)}%, consensus AF &ge; {(defaultConsFreq * 100).toFixed(0)}%).
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-3 shrink-0">
          <div className="flex items-center gap-2 bg-slate-50 border border-slate-200 rounded-xl px-3 py-1.5 text-xs text-slate-700">
            <SlidersHorizontal className="w-3.5 h-3.5 text-red-600" />
            <span>Min Breadth %:</span>
            <input
              type="number"
              min={0}
              max={100}
              step={5}
              value={minBreadthCutoff}
              onChange={(e) => setMinBreadthCutoff(Math.max(0, Math.min(100, Number(e.target.value) || 0)))}
              className="w-14 bg-white border border-slate-300 rounded px-1.5 py-0.5 text-xs font-mono text-slate-900"
            />
          </div>

          <div className="flex items-center gap-2 bg-slate-50 border border-slate-200 rounded-xl px-3 py-1.5 text-xs text-slate-700">
            <span>Clonal AF %:</span>
            <input
              type="number"
              min={10}
              max={100}
              step={5}
              value={clonalAfPct}
              onChange={(e) => setClonalAfPct(Math.max(10, Math.min(100, Number(e.target.value) || 50)))}
              className="w-14 bg-white border border-slate-300 rounded px-1.5 py-0.5 text-xs font-mono text-slate-900"
            />
          </div>

          <div className="flex bg-slate-100 border border-slate-200 rounded-lg p-0.5 text-xs">
            <button
              onClick={() => setVariantScope('ALL')}
              className={`px-2.5 py-1 rounded-md font-medium transition-colors ${
                variantScope === 'ALL' ? 'bg-red-600 text-white font-semibold' : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              All Variants
            </button>
            <button
              onClick={() => setVariantScope('CONSENSUS')}
              className={`px-2.5 py-1 rounded-md font-medium transition-colors ${
                variantScope === 'CONSENSUS' ? 'bg-red-600 text-white font-semibold' : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              Consensus Only
            </button>
          </div>

          {refs && refs.length > 1 && setSelectedRefName && (
            <select
              value={selectedRefName}
              onChange={(e) => setSelectedRefName(e.target.value)}
              className="bg-slate-50 border border-slate-200 text-slate-900 text-xs rounded-xl px-3 py-1.5 focus:outline-none cursor-pointer"
            >
              {refs.map((r) => (
                <option key={r.name} value={r.name}>
                  {r.name}
                </option>
              ))}
            </select>
          )}
        </div>
      </div>

      {passingSamples.length === 0 ? (
        <div className="p-8 text-center text-slate-500 border border-slate-200 rounded-xl bg-slate-50 flex flex-col items-center gap-2">
          <AlertCircle className="w-8 h-8 text-amber-500" />
          <p className="font-semibold text-slate-900">
            No samples have sufficient coverage (&ge; {minBreadthCutoff}%) for reference: {selectedRefName || 'all'}
          </p>
          <p className="text-xs text-slate-500">Lower the Min Breadth filter above or select another reference strain.</p>
        </div>
      ) : sortedMutations.length === 0 ? (
        <div className="p-8 text-center text-slate-500 border border-slate-200 rounded-xl bg-slate-50">
          <p>
            No {variantScope === 'CONSENSUS' ? 'consensus ' : ''}variants identified across the {passingSamples.length} passing samples for {selectedRefName || 'all'}.
          </p>
        </div>
      ) : (
        <div className="overflow-x-auto border border-slate-200 rounded-xl bg-white max-w-full">
          <table className="border-collapse w-full min-w-[700px] text-xs">
            <thead>
              <tr className="border-b border-slate-200 bg-slate-100">
                <th
                  className="p-3 text-left font-semibold text-slate-700 border-r border-slate-200 w-48 shrink-0 cursor-pointer hover:text-slate-900 select-none"
                  onClick={() => setSortAsc(!sortAsc)}
                >
                  <div className="flex items-center gap-1">
                    Sample ID <ArrowUpDown className="w-3 h-3 text-red-600" />
                  </div>
                </th>
                {sortedMutations.map((mut) => {
                  const primaryCons = mut.consequences[0] || null;
                  return (
                    <th key={mut.id} className="p-2.5 text-center font-mono font-medium text-slate-700 border-r border-slate-200 min-w-20">
                      <div className="flex flex-col items-center">
                        <span className="text-[10px] text-indigo-600 font-bold">{mut.gene}</span>
                        <span className="text-slate-900 font-semibold mt-0.5">{mut.ref}{mut.position}{mut.alt}</span>
                        {primaryCons && (
                          <span className={`text-[9px] mt-0.5 ${
                            primaryCons.kind === 'unavailable' || primaryCons.kind === 'synonymous'
                              ? 'text-slate-400 font-normal'
                              : primaryCons.kind === 'frameshift' || primaryCons.kind === 'indel'
                              ? 'text-amber-700 font-semibold'
                              : 'text-red-600 font-bold'
                          }`}>
                            {mut.consequences.length > 1
                              ? mut.consequences.map((c) => c.label).join(' / ')
                              : primaryCons.label}
                          </span>
                        )}
                      </div>
                    </th>
                  );
                })}
              </tr>
            </thead>

            <tbody>
              {sortedPassingSamples.map((sample) => (
                <tr key={sample.sample_id} className="border-b border-slate-100 bg-white hover:bg-slate-50 transition-colors">
                  <td className="p-3 font-semibold text-slate-900 border-r border-slate-200 truncate max-w-48">
                    {sample.sample_id}
                  </td>

                  {sortedMutations.map((mut) => {
                    const cell = variantLookup.get(sample.sample_id)?.get(mut.id) || null;
                    const hasMut = cell !== null;
                    const isHighFreq = hasMut && cell.frequency >= clonalAf;
                    const { isLowCov } = hasMut
                      ? { isLowCov: false }
                      : isVariantPositionLowCoverage(sample, selectedRefName, mut.position);

                    return (
                      <td
                        key={mut.id}
                        className="p-1 border-r border-slate-100 text-center"
                        onMouseEnter={() => setHoveredCell({ sampleId: sample.sample_id, mut, cell, isLowCov })}
                        onMouseLeave={() => setHoveredCell(null)}
                      >
                        <div
                          className={`mx-auto w-8 h-8 rounded-lg flex items-center justify-center transition-all duration-150 ${
                            hasMut
                              ? isHighFreq
                                ? 'bg-red-600 text-white font-bold scale-95 border border-red-700'
                                : 'bg-red-100 text-red-700 border border-red-300 font-semibold scale-95'
                              : isLowCov
                              ? 'bg-slate-200 border border-slate-300 text-slate-600 scale-90'
                              : 'bg-white border border-slate-200 hover:border-slate-400'
                          }`}
                        >
                          {hasMut ? (
                            <span className="text-[10px] font-mono">{Math.round(cell.frequency * 100)}</span>
                          ) : isLowCov ? (
                            <span className="text-[9px] font-mono font-bold">N</span>
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

      {passingSamples.length > 0 && sortedMutations.length > 0 && (
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4 items-center">
          <div className="p-4 bg-slate-50 border border-slate-200 rounded-xl min-h-20">
            {hoveredCell ? (
              <div className="flex flex-col gap-1">
                <span className="text-[10px] text-slate-500 font-mono">{hoveredCell.sampleId}</span>
                <div className="flex flex-wrap items-center gap-2 mt-0.5">
                  <span className="text-sm font-bold text-slate-900 font-mono">
                    g.{hoveredCell.mut.position}{hoveredCell.mut.ref}&gt;{hoveredCell.mut.alt}
                  </span>
                  <span className="text-[10px] bg-indigo-50 text-indigo-700 border border-indigo-200 font-bold px-1.5 py-0.5 rounded">
                    {hoveredCell.mut.gene}
                  </span>
                  {hoveredCell.mut.consequences.map((c, cIdx) => (
                    <span
                      key={cIdx}
                      className={`text-[10px] font-bold px-1.5 py-0.5 rounded border ${
                        c.kind === 'unavailable' || c.kind === 'synonymous'
                          ? 'bg-slate-100 text-slate-500 border-slate-200'
                          : c.kind === 'frameshift' || c.kind === 'indel'
                          ? 'bg-amber-50 text-amber-700 border-amber-200'
                          : 'bg-red-50 text-red-700 border-red-200'
                      }`}
                    >
                      {hoveredCell.mut.consequences.length > 1 ? `${c.gene}: ` : ''}
                      {c.label}
                    </span>
                  ))}
                </div>
                {hoveredCell.cell ? (
                  <div className="flex items-center gap-4 text-xs mt-1 text-slate-700">
                    <span>Allele Freq: <strong className="text-red-600 font-mono">{(hoveredCell.cell.frequency * 100).toFixed(1)}%</strong></span>
                    <span>Depth: <strong className="text-slate-900 font-mono">{hoveredCell.cell.depth}x</strong></span>
                  </div>
                ) : hoveredCell.isLowCov ? (
                  <div className="text-xs mt-1 text-red-600 font-medium">
                    Low Coverage / Masked &apos;N&apos; at ref position {hoveredCell.mut.position} (&lt;{minCoverageThreshold}x depth)
                  </div>
                ) : (
                  <div className="text-xs mt-1 text-emerald-700 font-medium">
                    Wild Type Reference Allele ({hoveredCell.mut.ref}) &ge;{minCoverageThreshold}x covered
                  </div>
                )}
              </div>
            ) : (
              <div className="flex items-center gap-2 text-slate-500 text-xs h-full py-2 justify-center">
                <HelpCircle className="w-4 h-4 text-red-600" />
                <span>Hover over grid cells to inspect indel-safe reference coverage and allele frequency</span>
              </div>
            )}
          </div>

          <div className="flex flex-wrap gap-4 justify-start md:justify-end text-xs text-slate-600">
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-red-600 border border-red-700" />
              <span>Clonal (&ge;{clonalAfPct}%)</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-red-100 border border-red-300" />
              <span>Subclonal (&lt;{clonalAfPct}%)</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-white border border-slate-300" />
              <span>Reference Match</span>
            </div>
            <div className="flex items-center gap-2">
              <div className="w-4 h-4 rounded bg-slate-200 border border-slate-300 flex items-center justify-center text-[9px] font-mono text-slate-700 font-bold">
                N
              </div>
              <span>Low Coverage (&lt;{minCoverageThreshold}x)</span>
            </div>
          </div>
        </div>
      )}

      {excludedSamples.length > 0 && (
        <div className="border border-slate-200 bg-slate-50 rounded-xl p-4 flex flex-col gap-2">
          <span className="text-xs font-semibold text-red-700 tracking-wider block">
            Excluded Low-Coverage Samples (&lt;{minBreadthCutoff}% Breadth at &ge;{minCoverageThreshold}x)
          </span>
          <p className="text-[11px] text-slate-600">
            Samples below the breadth threshold are excluded from matrix comparisons:
          </p>
          <div className="flex flex-wrap gap-2 mt-1">
            {excludedSamples.map(({ sample, breadth }) => (
              <span key={sample.sample_id} className="text-[11px] bg-white border border-slate-200 text-slate-700 font-mono px-2.5 py-1 rounded-lg">
                {sample.sample_id}: <strong className="text-red-600">{breadth.toFixed(1)}% cov</strong>
              </span>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
