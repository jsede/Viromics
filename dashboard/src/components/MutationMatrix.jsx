import React, { useMemo, useState } from 'react';
import { AlertCircle, ArrowUpDown, SlidersHorizontal } from 'lucide-react';
import { deduplicateConsequences, formatConsequences, formatDepth, getEvalBreadth } from '../utils/reportSchema';
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
            consequences: deduplicateConsequences(formatConsequences(v)),
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
    <div className="p-5 glass-panel flex flex-col gap-4">
      <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-3">
        <div>
          <h3 className="text-[16px] font-semibold text-[#111827]">Mutational Sharing Matrix</h3>
          <p className="text-[12px] text-[#6B7280] mt-0.5 tabular-nums">
            Reference-coordinate sharing with indel-safe dropout lookup (AF &ge; {(defaultMinFreq * 100).toFixed(0)}%, consensus AF &ge; {(defaultConsFreq * 100).toFixed(0)}%)
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-2.5 shrink-0">
          <div className="flex items-center gap-1.5 bg-white border border-[#CBD2D9] rounded-[4px] px-2.5 py-1 text-xs text-[#4B5563]">
            <SlidersHorizontal className="w-3.5 h-3.5 text-[#0F5C6E]" />
            <span>Min Breadth %:</span>
            <input
              type="number"
              min={0}
              max={100}
              step={5}
              value={minBreadthCutoff}
              onChange={(e) => setMinBreadthCutoff(Math.max(0, Math.min(100, Number(e.target.value) || 0)))}
              className="w-14 bg-[#F8F9FA] border border-[#CBD2D9] rounded-[3px] px-1.5 py-0.5 text-xs font-mono text-[#111827] tabular-nums"
            />
          </div>

          <div className="flex items-center gap-1.5 bg-white border border-[#CBD2D9] rounded-[4px] px-2.5 py-1 text-xs text-[#4B5563]">
            <span>Clonal AF %:</span>
            <input
              type="number"
              min={10}
              max={100}
              step={5}
              value={clonalAfPct}
              onChange={(e) => setClonalAfPct(Math.max(10, Math.min(100, Number(e.target.value) || 50)))}
              className="w-14 bg-[#F8F9FA] border border-[#CBD2D9] rounded-[3px] px-1.5 py-0.5 text-xs font-mono text-[#111827] tabular-nums"
            />
          </div>

          <div className="flex bg-[#F1F3F5] border border-[#E2E6EA] rounded-[4px] p-0.5 text-xs">
            <button
              onClick={() => setVariantScope('ALL')}
              className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                variantScope === 'ALL' ? 'bg-[#0F5C6E] text-white font-semibold' : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              All Variants
            </button>
            <button
              onClick={() => setVariantScope('CONSENSUS')}
              className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                variantScope === 'CONSENSUS' ? 'bg-[#0F5C6E] text-white font-semibold' : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              Consensus Only
            </button>
          </div>

          {refs && refs.length > 1 && setSelectedRefName && (
            <select
              value={selectedRefName}
              onChange={(e) => setSelectedRefName(e.target.value)}
              className="bg-white border border-[#CBD2D9] text-[#111827] text-xs font-mono rounded-[4px] px-2.5 py-1.5 focus:outline-none focus:border-[#0F5C6E] cursor-pointer"
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

      {passingSamples.length > 0 && sortedMutations.length > 0 && (
        <div className="flex flex-wrap items-center gap-x-5 gap-y-1.5 text-xs text-[#4B5563] border-t border-[#E2E6EA] pt-2.5 tabular-nums">
          <div className="flex items-center gap-2 whitespace-nowrap">
            <div className="w-3.5 h-3.5 rounded-[3px] bg-[#0F5C6E] border border-[#0B4654]" />
            <span>Clonal (&ge;{clonalAfPct}%)</span>
          </div>
          <div className="flex items-center gap-2 whitespace-nowrap">
            <div className="w-3.5 h-3.5 rounded-[3px] bg-[#D4E8ED] border border-dashed border-[#0F5C6E]" />
            <span>Subclonal (&lt;{clonalAfPct}%)</span>
          </div>
          <div className="flex items-center gap-2 whitespace-nowrap">
            <div className="w-3.5 h-3.5 rounded-[3px] bg-white border border-[#CBD2D9]" />
            <span>Reference Match</span>
          </div>
          <div className="flex items-center gap-2 whitespace-nowrap">
            <div className="w-3.5 h-3.5 rounded-[3px] bg-[#E5E7EB] border border-[#CBD2D9] flex items-center justify-center text-[9px] font-mono text-[#4B5563] font-bold">
              N
            </div>
            <span>Low Coverage (&lt;{minCoverageThreshold}x)</span>
          </div>
        </div>
      )}

      {passingSamples.length === 0 ? (
        <div className="p-8 text-center text-[#6B7280] border border-[#E2E6EA] rounded-[6px] bg-[#F8F9FA] flex flex-col items-center gap-2">
          <AlertCircle className="w-6 h-6 text-[#8C5308]" />
          <p className="font-semibold text-[#111827]">
            No samples have sufficient coverage (&ge; {minBreadthCutoff}%) for reference: {selectedRefName || 'all'}
          </p>
          <p className="text-xs text-[#6B7280]">Lower the Min Breadth filter above or select another reference strain.</p>
        </div>
      ) : sortedMutations.length === 0 ? (
        <div className="p-8 text-center text-[#6B7280] border border-[#E2E6EA] rounded-[6px] bg-[#F8F9FA] text-[13px]">
          <p>
            No {variantScope === 'CONSENSUS' ? 'consensus ' : ''}variants identified across the {passingSamples.length} passing samples for {selectedRefName || 'all'}.
          </p>
        </div>
      ) : (
        <div className="overflow-x-auto border border-[#E2E6EA] rounded-[6px] bg-white max-w-full">
          <table className="border-collapse w-full min-w-[700px] text-xs">
            <thead>
              <tr className="border-b border-[#E2E6EA] bg-[#F1F3F5]">
                <th
                  className="p-3 text-left font-semibold text-[#4B5563] border-r border-[#E2E6EA] w-48 shrink-0 cursor-pointer hover:text-[#111827] select-none text-[11px] uppercase tracking-[0.04em]"
                  onClick={() => setSortAsc(!sortAsc)}
                >
                  <div className="flex items-center gap-1">
                    Sample ID <ArrowUpDown className="w-3 h-3 text-[#0F5C6E]" />
                  </div>
                </th>
                {sortedMutations.map((mut) => {
                  const primaryCons = mut.consequences[0] || null;
                  return (
                    <th key={mut.id} className="p-2 text-center font-mono font-medium text-[#4B5563] border-r border-[#E2E6EA] last:border-r-0 min-w-20">
                      <div className="flex flex-col items-center">
                        <span className="text-[10px] text-[#6B7280] font-semibold">{mut.gene}</span>
                        <span className="text-[#111827] font-semibold mt-0.5 tabular-nums">
                          {mut.ref}{mut.position}{mut.alt}
                        </span>
                        {primaryCons && (
                          <span
                            className={`text-[10px] mt-0.5 ${
                              primaryCons.kind === 'unavailable' || primaryCons.kind === 'synonymous'
                                ? 'text-[#6B7280] font-normal'
                                : primaryCons.kind === 'frameshift' || primaryCons.kind === 'indel'
                                ? 'text-[#9A5B06] font-semibold'
                                : 'text-[#B43B2E] font-semibold'
                            }`}
                          >
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
                <tr key={sample.sample_id} className="border-b border-[#E2E6EA] last:border-none bg-white hover:bg-[#F8F9FA] transition-colors">
                  <td className="px-3 py-2 font-semibold text-[#111827] border-r border-[#E2E6EA] truncate max-w-48">
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
                        className="p-1 border-r border-[#E2E6EA] last:border-r-0 text-center"
                        onMouseEnter={() => setHoveredCell({ sampleId: sample.sample_id, mut, cell, isLowCov })}
                        onMouseLeave={() => setHoveredCell(null)}
                      >
                        <div
                          className={`mx-auto w-8 h-8 rounded-[4px] flex items-center justify-center transition-colors font-mono tabular-nums ${
                            hasMut
                              ? isHighFreq
                                ? 'bg-[#0F5C6E] text-white font-semibold border border-[#0B4654]'
                                : 'bg-[#D4E8ED] text-[#0B4654] border border-dashed border-[#0F5C6E] font-semibold'
                              : isLowCov
                              ? 'bg-[#E5E7EB] border border-[#CBD2D9] text-[#4B5563]'
                              : 'bg-white border border-[#E2E6EA] hover:border-[#0F5C6E]'
                          }`}
                        >
                          {hasMut ? (
                            <span className="text-[10px]">{Math.round(cell.frequency * 100)}</span>
                          ) : isLowCov ? (
                            <span className="text-[10px] font-bold">N</span>
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
        <div className="px-3.5 py-2.5 bg-[#F8F9FA] border border-[#E2E6EA] rounded-[4px] min-h-[42px] flex flex-wrap items-center justify-between gap-2 text-xs">
          {hoveredCell ? (
            <>
              <div className="flex flex-wrap items-center gap-2">
                <span className="font-semibold text-[#111827]">{hoveredCell.sampleId}</span>
                <span className="text-[#CBD2D9]">·</span>
                <span className="font-mono font-semibold text-[#111827]">
                  g.{hoveredCell.mut.position}{hoveredCell.mut.ref}&gt;{hoveredCell.mut.alt}
                </span>
                <span className="text-[10px] bg-white text-[#4B5563] border border-[#E2E6EA] font-mono px-1.5 py-0.5 rounded-[4px]">
                  {hoveredCell.mut.gene}
                </span>
                {hoveredCell.mut.consequences.map((c, cIdx) => (
                  <span
                    key={cIdx}
                    className={`text-[10px] font-mono px-1.5 py-0.5 rounded-[4px] border ${
                      c.kind === 'unavailable' || c.kind === 'synonymous'
                        ? 'bg-[#F1F3F5] text-[#5C6773] border-[#E2E6EA]'
                        : c.kind === 'frameshift' || c.kind === 'indel'
                        ? 'bg-[#FDF6E7] text-[#9A5B06] border-[#EBD29B] font-semibold'
                        : 'bg-[#FCEEEB] text-[#B43B2E] border-[#EBBDB8] font-semibold'
                    }`}
                  >
                    {hoveredCell.mut.consequences.length > 1 ? `${c.gene}: ` : ''}
                    {c.label}
                  </span>
                ))}
              </div>
              {hoveredCell.cell ? (
                <div className="flex items-center gap-3 font-mono tabular-nums text-[#4B5563]">
                  <span>
                    Allele Freq: <strong className="text-[#0F5C6E]">{(hoveredCell.cell.frequency * 100).toFixed(1)}%</strong>
                  </span>
                  <span>
                    Depth: <strong className="text-[#111827]">{formatDepth(hoveredCell.cell.depth, { decimals: 0 })}</strong>
                  </span>
                </div>
              ) : hoveredCell.isLowCov ? (
                <span className="text-[#A62C22] font-medium">
                  Masked &apos;N&apos; at ref position {hoveredCell.mut.position} (&lt;{minCoverageThreshold}x depth)
                </span>
              ) : (
                <span className="text-[#16654A] font-medium">
                  Wild-type reference allele ({hoveredCell.mut.ref}) &ge;{minCoverageThreshold}x covered
                </span>
              )}
            </>
          ) : (
            <span className="text-[#6B7280]">
              Hover over any matrix cell to inspect indel-safe reference coverage, depth, and allele frequency.
            </span>
          )}
        </div>
      )}

      {excludedSamples.length > 0 && (
        <div className="border border-[#E2E6EA] bg-[#F8F9FA] rounded-[4px] px-4 py-3 flex flex-col gap-1.5">
          <span className="text-xs font-semibold text-[#A62C22]">
            Excluded Low-Coverage Samples (&lt;{minBreadthCutoff}% Breadth at &ge;{minCoverageThreshold}x)
          </span>
          <div className="flex flex-wrap gap-2">
            {excludedSamples.map(({ sample, breadth }) => (
              <span key={sample.sample_id} className="text-[11px] bg-white border border-[#E2E6EA] text-[#111827] font-mono px-2.5 py-1 rounded-[4px] tabular-nums">
                {sample.sample_id}: <strong className="text-[#A62C22]">{breadth.toFixed(1)}% cov</strong>
              </span>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
