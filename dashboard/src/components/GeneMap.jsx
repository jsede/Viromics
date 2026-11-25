import React, { useMemo, useRef, useState } from 'react';
import { HelpCircle, Info, Download } from 'lucide-react';
import { formatConsequences } from '../utils/reportSchema';
import { serializeSvgElement, triggerBrowserDownload } from '../utils/exporters';

const FEATURE_PALETTE = [
  { fill: '#DC2626', stroke: '#B91C1C', badge: 'bg-red-100 text-red-700' },
  { fill: '#4F46E5', stroke: '#4338CA', badge: 'bg-indigo-100 text-indigo-700' },
  { fill: '#059669', stroke: '#047857', badge: 'bg-emerald-100 text-emerald-700' },
  { fill: '#7C3AED', stroke: '#6D28D9', badge: 'bg-violet-100 text-violet-700' },
  { fill: '#D97706', stroke: '#B45309', badge: 'bg-amber-100 text-amber-700' },
  { fill: '#0284C7', stroke: '#0369A1', badge: 'bg-sky-100 text-sky-700' },
];

function getFeatureDisplayName(gene) {
  const name = gene.name || 'Unknown';
  const prod = (gene.product || '').toLowerCase();
  if (name === 'ORF1ab' && prod.includes('orf1a polyprotein')) return 'ORF1a';
  return name;
}

function packFeatureTracks(genes) {
  if (!Array.isArray(genes) || genes.length === 0) return { packed: [], numTracks: 1 };

  // Group multi-part CDS (e.g. ORF1ab) so exon segments share the same visual track & connector
  const groupsMap = new Map();
  genes.forEach((g, idx) => {
    const key = g.grouping_key
      ? `${g.strand || '+'}::${g.grouping_key}`
      : `single_${idx}`;
    if (!groupsMap.has(key)) {
      groupsMap.set(key, {
        key,
        label: getFeatureDisplayName(g),
        strand: g.strand || '+',
        product: g.product || '',
        protein_id: g.protein_id || '',
        parts: [],
        minStart: g.start,
        maxEnd: g.end,
      });
    }
    const grp = groupsMap.get(key);
    grp.parts.push(g);
    grp.minStart = Math.min(grp.minStart, g.start);
    grp.maxEnd = Math.max(grp.maxEnd, g.end);
  });

  const groups = [...groupsMap.values()].sort((a, b) => a.minStart - b.minStart || b.maxEnd - a.maxEnd);
  const trackEnds = [];
  const packed = groups.map((grp, colorIdx) => {
    let assignedTrack = trackEnds.findIndex((endPos) => grp.minStart > endPos);
    if (assignedTrack === -1) {
      assignedTrack = trackEnds.length;
      trackEnds.push(grp.maxEnd);
    } else {
      trackEnds[assignedTrack] = grp.maxEnd;
    }
    return {
      ...grp,
      track: assignedTrack,
      color: FEATURE_PALETTE[colorIdx % FEATURE_PALETTE.length],
    };
  });

  return { packed, numTracks: Math.max(1, trackEnds.length) };
}

export default function GeneMap({ genes, variants, refLength, sampleId = 'sample', selectedRefName = 'contig' }) {
  const [hoveredGroup, setHoveredGroup] = useState(null);
  const [hoveredVariant, setHoveredVariant] = useState(null);
  const svgRef = useRef(null);

  const { packed, numTracks } = useMemo(() => packFeatureTracks(genes || []), [genes]);

  if (!genes || genes.length === 0) {
    return (
      <div className="p-6 glass-panel rounded-2xl text-center text-slate-500">
        <Info className="w-6 h-6 mx-auto mb-2 text-red-600" />
        <p>No annotation (GFF3) file loaded. Custom genome browser disabled.</p>
      </div>
    );
  }

  const safeRefLength = refLength && refLength > 0
    ? refLength
    : Math.max(...genes.map((g) => g.end), 1);

  const trackHeight = 30;
  const topPad = 24;
  const axisY = topPad + numTracks * trackHeight + 16;
  const svgHeight = axisY + 65;

  const handleExportSvg = () => {
    if (!svgRef.current) return;
    const svgText = serializeSvgElement(svgRef.current);
    triggerBrowserDownload(svgText, `viromics_${sampleId}_${selectedRefName}_genemap.svg`, 'image/svg+xml;charset=utf-8;');
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-5">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3">
        <div>
          <h3 className="text-lg font-semibold text-slate-900">Genome Feature Browser</h3>
          <p className="text-xs text-slate-500">
            Non-overlapping staggered track view of GFF3 CDS/gene features with mapped variant positions
          </p>
        </div>
        <button
          onClick={handleExportSvg}
          className="self-start sm:self-auto flex items-center gap-1.5 px-3 py-1.5 bg-white border border-slate-300 text-slate-700 text-xs font-semibold rounded-lg hover:bg-slate-50 transition-colors"
        >
          <Download className="w-3.5 h-3.5 text-red-600" />
          Export SVG
        </button>
      </div>

      <div className="relative w-full overflow-x-auto py-2 select-none">
        <svg ref={svgRef} viewBox={`0 0 1000 ${svgHeight}`} width="100%" className="min-w-[800px]">
          <line x1="20" y1={axisY} x2="980" y2={axisY} stroke="#CBD5E1" strokeWidth="5" strokeLinecap="round" />

          {[0, 0.25, 0.5, 0.75, 1].map((pct, i) => {
            const pos = Math.round(pct * safeRefLength);
            const x = 20 + pct * 960;
            return (
              <g key={i}>
                <line x1={x} y1={axisY} x2={x} y2={axisY + 6} stroke="#94A3B8" strokeWidth="1.5" />
                <text x={x} y={axisY + 52} fill="#64748B" fontSize="9.5" textAnchor="middle" fontFamily="Fira Code, monospace">
                  {pos.toLocaleString()} bp
                </text>
              </g>
            );
          })}

          {packed.map((grp) => {
            const y = topPad + grp.track * trackHeight;
            const h = 13;
            const isHovered = hoveredGroup && hoveredGroup.key === grp.key;
            const spanX = 20 + (grp.minStart / safeRefLength) * 960;
            const spanW = Math.max(8, ((grp.maxEnd - grp.minStart) / safeRefLength) * 960);

            return (
              <g
                key={grp.key}
                onMouseEnter={() => setHoveredGroup(grp)}
                onMouseLeave={() => setHoveredGroup(null)}
                className="cursor-pointer"
              >
                {grp.parts.length > 1 && (
                  <line
                    x1={spanX}
                    y1={y + h / 2}
                    x2={spanX + spanW}
                    y2={y + h / 2}
                    stroke={grp.color.stroke}
                    strokeWidth="1.5"
                    strokeDasharray="3 2"
                  />
                )}
                {grp.parts.map((part, pIdx) => {
                  const x = 20 + (part.start / safeRefLength) * 960;
                  const width = Math.max(6, ((part.end - part.start) / safeRefLength) * 960);
                  return (
                    <rect
                      key={`${grp.key}-p-${pIdx}`}
                      x={x}
                      y={y}
                      width={width}
                      height={h}
                      rx="3"
                      fill={grp.color.fill}
                      fillOpacity={isHovered ? 0.85 : 0.4}
                      stroke={grp.color.stroke}
                      strokeOpacity={isHovered ? 1 : 0.75}
                      strokeWidth="1.5"
                      className="gene-block"
                    />
                  );
                })}
                {spanW > 24 && (
                  <text
                    x={spanX + spanW / 2}
                    y={y - 5}
                    fill={isHovered ? '#DC2626' : '#334155'}
                    fontSize="9.5"
                    fontWeight={isHovered ? '700' : '600'}
                    textAnchor="middle"
                  >
                    {grp.label}
                  </text>
                )}
              </g>
            );
          })}

          {variants &&
            variants.map((v, idx) => {
              const x = 20 + (v.position / safeRefLength) * 960;
              const isHovered = hoveredVariant && hoveredVariant.position === v.position && hoveredVariant.alt === v.alt;
              return (
                <g
                  key={idx}
                  onMouseEnter={() => setHoveredVariant(v)}
                  onMouseLeave={() => setHoveredVariant(null)}
                  className="cursor-pointer"
                >
                  <line
                    x1={x}
                    y1={axisY}
                    x2={x}
                    y2={axisY + 28}
                    stroke={isHovered ? '#DC2626' : '#64748B'}
                    strokeWidth={isHovered ? 2.5 : 1.5}
                  />
                  <circle
                    cx={x}
                    cy={axisY + 28}
                    r={isHovered ? 5.5 : 3.5}
                    fill={isHovered ? '#DC2626' : '#64748B'}
                  />
                </g>
              );
            })}
        </svg>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
        <div className="p-4 bg-slate-50 border border-slate-200 rounded-xl min-h-24">
          {hoveredGroup ? (
            <div className="flex flex-col gap-1.5">
              <div className="flex items-center gap-2">
                <span className={`px-2 py-0.5 rounded text-[10px] font-bold ${hoveredGroup.color.badge}`}>
                  CDS ({hoveredGroup.parts.length} {hoveredGroup.parts.length === 1 ? 'interval' : 'intervals'})
                </span>
                <h4 className="text-sm font-bold text-slate-900">{hoveredGroup.label}</h4>
              </div>
              <p className="text-xs text-slate-700 font-medium">
                {hoveredGroup.product || 'No product annotation available'}
              </p>
              <div className="flex flex-wrap items-center gap-4 text-[10px] font-mono text-slate-500 mt-0.5">
                <span>Span: {hoveredGroup.minStart.toLocaleString()}–{hoveredGroup.maxEnd.toLocaleString()} bp</span>
                <span>Strand: ({hoveredGroup.strand})</span>
                {hoveredGroup.protein_id && <span>Protein: {hoveredGroup.protein_id}</span>}
              </div>
            </div>
          ) : (
            <div className="flex items-center gap-2 text-slate-500 text-xs h-full py-4 justify-center">
              <HelpCircle className="w-4 h-4 text-indigo-600" />
              <span>Hover over a feature block to inspect CDS intervals and product details</span>
            </div>
          )}
        </div>

        <div className="p-4 bg-slate-50 border border-slate-200 rounded-xl min-h-24">
          {hoveredVariant ? (
            <div className="flex flex-col gap-1.5">
              <div className="flex items-center gap-2">
                <span className="px-2 py-0.5 rounded text-[10px] font-bold bg-red-100 text-red-700">
                  VARIANT
                </span>
                <h4 className="text-sm font-bold text-slate-900 font-mono">
                  g.{hoveredVariant.position}{hoveredVariant.ref}&gt;{hoveredVariant.alt}
                </h4>
              </div>
              <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs mt-0.5">
                <span className="text-slate-500">Gene Region: <strong className="text-slate-900">{hoveredVariant.gene}</strong></span>
                <span className="text-slate-500">Class: <strong className="text-slate-900">{hoveredVariant.type}</strong></span>
                <span className="text-slate-500">Allele Freq: <strong className="text-red-600 font-mono">{(hoveredVariant.frequency * 100).toFixed(1)}%</strong></span>
                <span className="text-slate-500">Depth: <strong className="text-slate-900 font-mono">{hoveredVariant.depth}x</strong></span>
              </div>
              {(() => {
                const cons = formatConsequences(hoveredVariant);
                if (cons.length === 0) return null;
                return (
                  <div className="flex flex-wrap gap-1.5 mt-1">
                    {cons.map((c, idx) => (
                      <span
                        key={idx}
                        className={`text-[10px] px-2 py-0.5 rounded border font-mono ${
                          c.kind === 'unavailable' || c.kind === 'synonymous'
                            ? 'bg-slate-100 text-slate-500 border-slate-200'
                            : c.kind === 'frameshift' || c.kind === 'indel'
                            ? 'bg-amber-50 text-amber-700 border-amber-200 font-semibold'
                            : 'bg-indigo-50 text-indigo-700 border-indigo-200 font-bold'
                        }`}
                      >
                        {c.gene}: {c.label}
                      </span>
                    ))}
                  </div>
                );
              })()}
            </div>
          ) : (
            <div className="flex items-center gap-2 text-slate-500 text-xs h-full py-4 justify-center">
              <HelpCircle className="w-4 h-4 text-red-600" />
              <span>Hover over a variant marker pin to inspect structured consequences</span>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
