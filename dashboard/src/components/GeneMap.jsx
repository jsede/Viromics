import React, { useState } from 'react';
import { HelpCircle, Info } from 'lucide-react';

export default function GeneMap({ genes, variants, refLength }) {
  const [hoveredGene, setHoveredGene] = useState(null);
  const [hoveredVariant, setHoveredVariant] = useState(null);

  if (!genes || genes.length === 0) {
    return (
      <div className="p-6 glass-panel rounded-2xl text-center text-gray-400">
        <Info className="w-6 h-6 mx-auto mb-2 text-brand-teal animate-pulse" />
        <p>No annotation (GFF) file loaded. Custom genome browser disabled.</p>
      </div>
    );
  }

  // Color mapping based on gene names (so key genes like Spike, Nucleocapsid stand out!)
  const getGeneColor = (name) => {
    const n = name.toLowerCase();
    if (n === 'spike' || n === 's') return 'fill-brand-teal stroke-brand-teal';
    if (n === 'nucleocapsid' || n === 'n') return 'fill-brand-violet stroke-brand-violet';
    if (n.startsWith('orf1')) return 'fill-brand-indigo/60 stroke-brand-indigo/70';
    if (n === 'envelope' || n === 'e') return 'fill-brand-amber stroke-brand-amber';
    if (n === 'membrane' || n === 'm') return 'fill-brand-emerald stroke-brand-emerald';
    return 'fill-gray-600 stroke-gray-500';
  };

  const getGeneBgColor = (name) => {
    const n = name.toLowerCase();
    if (n === 'spike' || n === 's') return 'bg-brand-teal/20 text-brand-teal';
    if (n === 'nucleocapsid' || n === 'n') return 'bg-brand-violet/20 text-brand-violet';
    if (n.startsWith('orf1')) return 'bg-brand-indigo/20 text-brand-indigo';
    if (n === 'envelope' || n === 'e') return 'bg-brand-amber/20 text-brand-amber';
    if (n === 'membrane' || n === 'm') return 'bg-brand-emerald/20 text-brand-emerald';
    return 'bg-gray-600/20 text-gray-400';
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      {/* Header */}
      <div className="flex items-center justify-between">
        <div>
          <h3 className="text-lg font-semibold font-outfit text-white">Genome Feature Browser</h3>
          <p className="text-xs text-gray-400">
            Linear view of genes annotated from GFF with mapped variant positions (vertical markers)
          </p>
        </div>
      </div>

      {/* SVG Canvas for Linear Map */}
      <div className="relative w-full overflow-x-auto py-4 select-none">
        <svg viewBox={`0 0 1000 160`} width="100%" className="min-w-[800px]">
          {/* Main Genome Axis */}
          <line 
            x1="20" y1="80" x2="980" y2="80" 
            stroke="#CBD5E1" strokeWidth="5" strokeLinecap="round" 
          />
          
          {/* Genome coordinate ticks */}
          {[0, 0.25, 0.5, 0.75, 1].map((pct, i) => {
            const pos = Math.round(pct * refLength);
            const x = 20 + pct * 960;
            return (
              <g key={i}>
                <line x1={x} y1={80} x2={x} y2={86} stroke="#94A3B8" strokeWidth="1.5" />
                <text 
                  x={x} y={145} 
                  fill="#64748B" fontSize="9" textAnchor="middle" 
                  fontFamily="Fira Code"
                >
                  {pos.toLocaleString()} bp
                </text>
              </g>
            );
          })}

          {/* Gene Blocks */}
          {genes.map((gene, idx) => {
            // Scale positions to [20, 980]
            const startPct = gene.start / refLength;
            const endPct = gene.end / refLength;
            const x = 20 + startPct * 960;
            const width = Math.max(8, (endPct - startPct) * 960);
            
            // Draw genes alternated top/bottom or height offset to avoid overlapping names
            const isAlternate = idx % 2 === 0;
            const y = isAlternate ? 20 : 48;
            const h = 12;
            
            const isHovered = hoveredGene && hoveredGene.name === gene.name;

            return (
              <g 
                key={idx}
                onMouseEnter={() => setHoveredGene(gene)}
                onMouseLeave={() => setHoveredGene(null)}
                className="cursor-pointer transition-all duration-200"
              >
                <rect 
                  x={x} y={y} width={width} height={h} rx="3"
                  className={`gene-block ${getGeneColor(gene.name)} ${
                    isHovered ? 'fill-opacity-85 stroke-opacity-100' : 'fill-opacity-40 stroke-opacity-60'
                  } stroke-[1.5]`}
                />
                
                {/* Gene Name Label inside/above block */}
                {width > 25 && (
                  <text 
                    x={x + width/2} y={y - 8}
                    fill={isHovered ? '#DC2626' : '#475569'}
                    fontSize="9.5" fontWeight={isHovered ? '700' : '600'}
                    textAnchor="middle"
                    className="transition-colors"
                  >
                    {gene.name}
                  </text>
                )}
              </g>
            );
          })}

          {/* Variant Marker Pins along the axis (Pointing Downwards to prevent overlaps) */}
          {variants && variants.map((v, idx) => {
            const pct = v.position / refLength;
            const x = 20 + pct * 960;
            const y = 80;
            
            const isHovered = hoveredVariant && hoveredVariant.position === v.position;

            return (
              <g
                key={idx}
                onMouseEnter={() => setHoveredVariant(v)}
                onMouseLeave={() => setHoveredVariant(null)}
                className="cursor-pointer"
              >
                {/* Vertical Pin Line */}
                <line 
                  x1={x} y1={y} x2={x} y2={115} 
                  stroke={isHovered ? '#DC2626' : '#94A3B8'} 
                  strokeWidth={isHovered ? 2.5 : 1.5}
                />
                {/* Pin Head */}
                <circle 
                  cx={x} cy={115} 
                  r={isHovered ? 5.5 : 3.5} 
                  fill={isHovered ? '#DC2626' : '#94A3B8'}
                  className="transition-all duration-150"
                />
              </g>
            );
          })}
        </svg>
      </div>

      {/* Tooltip detail sections */}
      <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mt-2">
        {/* Gene Tooltip Detail */}
        <div className="p-4 bg-dark-bg/60 border border-dark-border rounded-xl min-h-24 transition-all">
          {hoveredGene ? (
            <div className="flex flex-col gap-1.5">
              <div className="flex items-center gap-2">
                <span className={`px-2 py-0.5 rounded text-[10px] font-bold ${getGeneBgColor(hoveredGene.name)}`}>
                  GENE
                </span>
                <h4 className="text-sm font-bold text-white font-outfit">{hoveredGene.name}</h4>
              </div>
              <p className="text-xs text-gray-300 font-medium">{hoveredGene.product || 'No product product information'}</p>
              <div className="flex items-center gap-4 text-[10px] font-mono text-gray-400 mt-1">
                <span>Start: {hoveredGene.start.toLocaleString()} bp</span>
                <span>End: {hoveredGene.end.toLocaleString()} bp</span>
                <span>Strand: ({hoveredGene.strand})</span>
              </div>
            </div>
          ) : (
            <div className="flex items-center gap-2 text-gray-400 text-xs h-full py-4 justify-center">
              <HelpCircle className="w-4 h-4 text-brand-indigo" />
              <span>Hover over a gene block to inspect details</span>
            </div>
          )}
        </div>

        {/* Variant Tooltip Detail */}
        <div className="p-4 bg-dark-bg/60 border border-dark-border rounded-xl min-h-24 transition-all">
          {hoveredVariant ? (
            <div className="flex flex-col gap-1.5">
              <div className="flex items-center gap-2">
                <span className="px-2 py-0.5 rounded text-[10px] font-bold bg-brand-rose/20 text-brand-rose">
                  VARIANT
                </span>
                <h4 className="text-sm font-bold text-white font-mono">
                  g.{hoveredVariant.position}{hoveredVariant.ref}&gt;{hoveredVariant.alt}
                </h4>
              </div>
              <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs mt-1">
                <span className="text-gray-400">Gene Region: <strong className="text-gray-200">{hoveredVariant.gene}</strong></span>
                <span className="text-gray-400">Variant Class: <strong className="text-gray-200">{hoveredVariant.type}</strong></span>
                <span className="text-gray-400">Allele Freq: <strong className="text-brand-teal font-mono">{(hoveredVariant.frequency * 100).toFixed(1)}%</strong></span>
                <span className="text-gray-400">Coverage Depth: <strong className="text-gray-200 font-mono">{hoveredVariant.depth}x</strong></span>
                {hoveredVariant.aa_change && (
                  <span className="text-gray-400 col-span-2">
                    Protein Change:{' '}
                    <strong className={
                      hoveredVariant.aa_change.includes('(syn)') || hoveredVariant.aa_change === 'AA annotation unavailable'
                        ? 'text-gray-400 font-medium'
                        : hoveredVariant.aa_change === 'Indel'
                        ? 'text-brand-amber font-medium'
                        : 'text-brand-indigo font-bold'
                    }>
                      {hoveredVariant.aa_change.includes('(syn)') 
                        ? `${hoveredVariant.aa_change.replace(' (syn)', '')} (Synonymous)` 
                        : hoveredVariant.aa_change === 'AA annotation unavailable'
                        ? `AA unavailable${hoveredVariant.status ? ` (${hoveredVariant.status})` : ''}`
                        : hoveredVariant.aa_change === 'Indel'
                        ? 'Indel'
                        : `${hoveredVariant.aa_change.split('/').map(a => a.startsWith('p.') ? a : `p.${a}`).join('/')} (Non-synonymous)`}
                    </strong>
                  </span>
                )}
              </div>
            </div>
          ) : (
            <div className="flex items-center gap-2 text-gray-400 text-xs h-full py-4 justify-center">
              <HelpCircle className="w-4 h-4 text-brand-rose" />
              <span>Hover over a marker pin (vertical bar) to inspect mutation details</span>
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
