import React, { useMemo, useRef, useState } from 'react';
import { Info, Download } from 'lucide-react';
import { deduplicateConsequences, formatConsequences } from '../utils/reportSchema';
import { serializeSvgElement, triggerBrowserDownload } from '../utils/exporters';
import { getVariantConsequenceStyle } from '../utils/chartUtils';

const FEATURE_PALETTE = [
  { fill: '#33658A', stroke: '#244A66', badge: 'bg-[#EAF3F5] text-[#0B4654] border-[#B9D6DD]' },
  { fill: '#4A9386', stroke: '#356E64', badge: 'bg-[#ECF7F2] text-[#16654A] border-[#B5DFCE]' },
  { fill: '#B38B3A', stroke: '#8A6929', badge: 'bg-[#FDF6E7] text-[#8C5308] border-[#EBD29B]' },
  { fill: '#8B5E83', stroke: '#674360', badge: 'bg-[#F1F3F5] text-[#111827] border-[#CBD2D9]' },
  { fill: '#6B8E5B', stroke: '#4D6940', badge: 'bg-[#ECF7F2] text-[#16654A] border-[#B5DFCE]' },
  { fill: '#B86B52', stroke: '#8D4E39', badge: 'bg-[#FCEEEB] text-[#A62C22] border-[#EBBDB8]' },
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
    const key = g.grouping_key ? `${g.strand || '+'}::${g.grouping_key}` : `single_${idx}`;
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
      <div className="p-6 glass-panel text-center text-[#6B7280] text-[13px]">
        <Info className="w-5 h-5 mx-auto mb-2 text-[#0F5C6E]" />
        <p>No annotation (GFF3) file loaded. Custom genome browser disabled.</p>
      </div>
    );
  }

  const safeRefLength = refLength && refLength > 0 ? refLength : Math.max(...genes.map((g) => g.end), 1);

  const leftPad = 36;
  const rightPad = 36;
  const plotWidth = 1000 - leftPad - rightPad; // 928
  const trackHeight = 30;
  const topPad = 24;
  const axisY = topPad + numTracks * trackHeight + 16;
  const svgHeight = axisY + 66;

  const handleExportSvg = () => {
    if (!svgRef.current) return;
    const svgText = serializeSvgElement(svgRef.current);
    triggerBrowserDownload(svgText, `viromics_${sampleId}_${selectedRefName}_genemap.svg`, 'image/svg+xml;charset=utf-8;');
  };

  return (
    <div className="p-5 glass-panel flex flex-col gap-4">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3">
        <div>
          <h3 className="text-[16px] font-semibold text-[#111827]">Genome Feature Browser</h3>
          <p className="text-[12px] text-[#6B7280] mt-0.5">
            Staggered GFF3 CDS/gene tracks with consequence-encoded variant lollipops
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-3">
          <div className="flex flex-wrap items-center gap-3 text-[11px] text-[#4B5563]">
            <span className="flex items-center gap-1.5">
              <span className="w-2.5 h-2.5 rounded-full bg-[#B43B2E]" />
              Missense / Nonsense
            </span>
            <span className="flex items-center gap-1.5">
              <span className="w-2.5 h-2.5 rounded-full bg-[#CC8400]" />
              Indel / Frameshift
            </span>
            <span className="flex items-center gap-1.5">
              <span className="w-2.5 h-2.5 rounded-full bg-[#5C6773]" />
              Synonymous
            </span>
            <span className="flex items-center gap-1.5">
              <span className="w-2.5 h-2.5 rounded-full bg-[#6B7280]" />
              Non-coding
            </span>
          </div>

          <button
            onClick={handleExportSvg}
            className="flex items-center gap-1.5 px-3 py-1.5 bg-white border border-[#CBD2D9] text-[#111827] text-xs font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] transition-colors whitespace-nowrap"
          >
            <Download className="w-3.5 h-3.5 text-[#0F5C6E]" />
            Export SVG
          </button>
        </div>
      </div>

      <div className="relative w-full overflow-x-auto py-1 select-none">
        <svg
          ref={svgRef}
          viewBox={`0 0 1000 ${svgHeight}`}
          width="100%"
          className="min-w-[780px]"
          style={{ fontFamily: 'Arial, Helvetica, sans-serif' }}
        >
          <line
            x1={leftPad}
            y1={axisY}
            x2={leftPad + plotWidth}
            y2={axisY}
            stroke="#CBD2D9"
            strokeWidth="4"
            strokeLinecap="round"
          />

          {[0, 0.25, 0.5, 0.75, 1].map((pct, i) => {
            const pos = Math.round(pct * safeRefLength);
            const x = leftPad + pct * plotWidth;
            const textAnchor = i === 0 ? 'start' : i === 4 ? 'end' : 'middle';
            return (
              <g key={i}>
                <line x1={x} y1={axisY} x2={x} y2={axisY + 6} stroke="#6B7280" strokeWidth="1.25" />
                <text
                  x={x}
                  y={axisY + 52}
                  fill="#4B5563"
                  fontSize="10"
                  textAnchor={textAnchor}
                  fontFamily="Arial, Helvetica, sans-serif"
                >
                  {pos.toLocaleString()} bp
                </text>
              </g>
            );
          })}

          {packed.map((grp) => {
            const y = topPad + grp.track * trackHeight;
            const h = 13;
            const isHovered = hoveredGroup && hoveredGroup.key === grp.key;
            const spanX = leftPad + (grp.minStart / safeRefLength) * plotWidth;
            const spanW = Math.max(8, ((grp.maxEnd - grp.minStart) / safeRefLength) * plotWidth);

            return (
              <g
                key={grp.key}
                onMouseEnter={() => {
                  setHoveredVariant(null);
                  setHoveredGroup(grp);
                }}
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
                  const x = leftPad + (part.start / safeRefLength) * plotWidth;
                  const width = Math.max(6, ((part.end - part.start) / safeRefLength) * plotWidth);
                  return (
                    <rect
                      key={`${grp.key}-p-${pIdx}`}
                      x={x}
                      y={y}
                      width={width}
                      height={h}
                      rx="3"
                      fill={grp.color.fill}
                      fillOpacity={isHovered ? 0.85 : 0.42}
                      stroke={grp.color.stroke}
                      strokeOpacity={isHovered ? 1 : 0.8}
                      strokeWidth="1.25"
                      className="gene-block"
                    />
                  );
                })}
                {spanW > 24 && (
                  <text
                    x={spanX + spanW / 2}
                    y={y - 5}
                    fill={isHovered ? '#0F5C6E' : '#111827'}
                    fontSize="10"
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
              const x = leftPad + (v.position / safeRefLength) * plotWidth;
              const isHovered = hoveredVariant && hoveredVariant.position === v.position && hoveredVariant.alt === v.alt;
              const { color } = getVariantConsequenceStyle(v);
              return (
                <g
                  key={idx}
                  onMouseEnter={() => {
                    setHoveredGroup(null);
                    setHoveredVariant(v);
                  }}
                  onMouseLeave={() => setHoveredVariant(null)}
                  className="cursor-pointer"
                >
                  <line
                    x1={x}
                    y1={axisY}
                    x2={x}
                    y2={axisY + 26}
                    stroke={isHovered ? '#0F5C6E' : color}
                    strokeWidth={isHovered ? 2.25 : 1.5}
                  />
                  <circle
                    cx={x}
                    cy={axisY + 26}
                    r={isHovered ? 5 : 3.5}
                    fill={color}
                    stroke={isHovered ? '#0F5C6E' : '#FFFFFF'}
                    strokeWidth={isHovered ? 2 : 0.75}
                  />
                </g>
              );
            })}
        </svg>
      </div>

      <div className="px-3.5 py-2.5 bg-[#F8F9FA] border border-[#E2E6EA] rounded-[4px] min-h-[44px] flex flex-wrap items-center justify-between gap-2 text-xs">
        {hoveredVariant ? (
          <>
            <div className="flex flex-wrap items-center gap-2.5">
              <span className="px-2 py-0.5 rounded-[4px] text-[10px] font-semibold bg-[#EAF3F5] text-[#0B4654] border border-[#B9D6DD]">
                VARIANT
              </span>
              <span className="font-mono font-semibold text-[#111827]">
                g.{hoveredVariant.position}{hoveredVariant.ref}&gt;{hoveredVariant.alt}
              </span>
              <span className="text-[#4B5563]">
                Gene: <strong className="text-[#111827]">{hoveredVariant.gene}</strong>
              </span>
              <span className="text-[#4B5563] font-mono tabular-nums">
                AF: <strong className="text-[#0F5C6E]">{(hoveredVariant.frequency * 100).toFixed(1)}%</strong> · Depth:{' '}
                <strong className="text-[#111827]">{hoveredVariant.depth}x</strong>
              </span>
            </div>
            <div className="flex flex-wrap gap-1.5">
              {deduplicateConsequences(formatConsequences(hoveredVariant)).map((c, idx) => (
                <span
                  key={idx}
                  title={c.tooltip}
                  className={`text-[10px] px-2 py-0.5 rounded-[4px] border font-mono ${
                    c.kind === 'unavailable'
                      ? 'bg-[#F1F3F5] text-[#6B7280] border-[#E2E6EA] italic'
                      : c.kind === 'synonymous'
                      ? 'bg-[#F1F3F5] text-[#4B5563] border-[#E2E6EA]'
                      : c.kind === 'frameshift' || c.kind === 'indel'
                      ? 'bg-[#FDF6E7] text-[#9A5B06] border-[#EBD29B] font-semibold'
                      : 'bg-[#FCEEEB] text-[#B43B2E] border-[#EBBDB8] font-semibold'
                  }`}
                >
                  {c.gene}: {c.label}
                </span>
              ))}
            </div>
          </>
        ) : hoveredGroup ? (
          <>
            <div className="flex flex-wrap items-center gap-2.5">
              <span className={`px-2 py-0.5 rounded-[4px] border text-[10px] font-semibold ${hoveredGroup.color.badge}`}>
                CDS ({hoveredGroup.parts.length} {hoveredGroup.parts.length === 1 ? 'interval' : 'intervals'})
              </span>
              <strong className="text-[#111827] font-semibold">{hoveredGroup.label}</strong>
              <span className="text-[#4B5563]">{hoveredGroup.product || 'No product annotation available'}</span>
            </div>
            <div className="flex flex-wrap items-center gap-3 text-[11px] font-mono text-[#6B7280] tabular-nums">
              <span>
                Span: {hoveredGroup.minStart.toLocaleString()}–{hoveredGroup.maxEnd.toLocaleString()} bp
              </span>
              <span>Strand: ({hoveredGroup.strand})</span>
              {hoveredGroup.protein_id && <span>Protein: {hoveredGroup.protein_id}</span>}
            </div>
          </>
        ) : (
          <span className="text-[#6B7280]">
            Hover over a CDS track block or variant lollipop marker above to inspect coordinates, product annotation, and protein consequences.
          </span>
        )}
      </div>
    </div>
  );
}
