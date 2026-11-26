import React, { useMemo, useRef, useState } from 'react';
import { GitBranch, AlertCircle, Info, Download } from 'lucide-react';
import { buildNeighborJoiningTree, buildUPGMATree, filterTreeSamples } from '../utils/phylogenetics';
import {
  formatCohortMultiFasta,
  formatDistanceMatrixCsv,
  serializeSvgElement,
  treeToNewick,
  triggerBrowserDownload,
} from '../utils/exporters';

function computeTreeLayout(rootNode, numLeaves, maxLabelChars, leftMargin = 36, treeWidth = 540) {
  if (!rootNode) return { coords: {}, maxDepth: 0, svgWidth: 800, svgHeight: 260, leftMargin, treeWidth };

  const computeMaxDepth = (n) => {
    if (!n.children || n.children.length === 0) return n.depthFromRoot || 0;
    return Math.max(...n.children.map(computeMaxDepth));
  };
  const maxDepth = computeMaxDepth(rootNode);
  const labelMargin = Math.max(150, Math.min(360, Math.round(maxLabelChars * 7.2 + 30)));
  const svgWidth = leftMargin + treeWidth + labelMargin;
  const svgHeight = Math.max(260, 44 + numLeaves * 44 + 44);

  const coords = {};
  let leafIdx = 0;

  const assign = (n) => {
    const x =
      maxDepth > 0
        ? leftMargin + (n.depthFromRoot / maxDepth) * treeWidth
        : n.isLeaf
        ? leftMargin + Math.round(treeWidth * 0.25)
        : leftMargin;

    if (n.isLeaf || !n.children || n.children.length === 0) {
      const y = 36 + leafIdx * 44;
      leafIdx++;
      coords[n.id] = { x, y, node: n };
      return y;
    }

    const childYs = n.children.map(assign);
    const y = childYs.reduce((a, b) => a + b, 0) / childYs.length;
    coords[n.id] = { x, y, node: n };
    return y;
  };

  assign(rootNode);
  return { coords, maxDepth, svgWidth, svgHeight, leftMargin, treeWidth };
}

function collectEdges(node, coords, out = []) {
  if (!node || !node.children || node.children.length === 0) return out;
  const parent = coords[node.id];
  node.children.forEach((child) => {
    const ch = coords[child.id];
    if (parent && ch) {
      out.push({
        key: `${node.id}->${child.id}`,
        d: `M ${parent.x} ${parent.y} L ${parent.x} ${ch.y} L ${ch.x} ${ch.y}`,
        branchLength: child.branchLength || 0,
      });
    }
    collectEdges(child, coords, out);
  });
  return out;
}

export default function PhylogeneticTree({
  samples,
  selectedRefName,
  setSelectedRefName,
  refs = [],
  reportData = null,
  minCoverageThreshold = 10,
}) {
  const [algorithm, setAlgorithm] = useState('NJ'); // 'NJ' (midpoint-rooted default) | 'UPGMA'
  const [minBreadth, setMinBreadth] = useState(80);
  const [minComparedSites, setMinComparedSites] = useState(1);
  const [hoveredNode, setHoveredNode] = useState(null);
  const svgRef = useRef(null);

  const {
    includedSamples,
    excludedLowBreadth,
    excludedOverlap,
    distMatrix,
    snpMatrix,
    sitesMatrix,
  } = useMemo(
    () =>
      filterTreeSamples(samples, selectedRefName, reportData, {
        minBreadth,
        minComparedSites,
      }),
    [samples, selectedRefName, reportData, minBreadth, minComparedSites]
  );

  const treeData = useMemo(() => {
    if (includedSamples.length < 2) return null;
    const ids = includedSamples.map((s) => s.sample_id);
    return algorithm === 'NJ'
      ? buildNeighborJoiningTree(ids, distMatrix)
      : buildUPGMATree(ids, distMatrix);
  }, [includedSamples, distMatrix, algorithm]);

  const maxLabelChars = useMemo(
    () => Math.max(12, ...includedSamples.map((s) => s.sample_id.length)),
    [includedSamples]
  );

  const layout = useMemo(
    () => computeTreeLayout(treeData, includedSamples.length, maxLabelChars),
    [treeData, includedSamples.length, maxLabelChars]
  );

  const treeEdges = useMemo(
    () => (treeData ? collectEdges(treeData, layout.coords) : []),
    [treeData, layout.coords]
  );

  if (!samples || samples.length < 2) {
    return (
      <div className="p-8 text-center text-[#6B7280] glass-panel flex flex-col items-center gap-3 text-[13px]">
        <AlertCircle className="w-6 h-6 text-[#8C5308]" />
        <p>Phylogenetic clustering requires at least 2 samples in the run.</p>
      </div>
    );
  }

  const handleDownloadNewick = () => {
    const nwk = treeToNewick(treeData);
    if (!nwk) return;
    triggerBrowserDownload(nwk, `viromics_${selectedRefName || 'tree'}_${algorithm.toLowerCase()}.nwk`, 'text/plain;charset=utf-8;');
  };

  const handleDownloadMatrixCsv = () => {
    const ids = includedSamples.map((s) => s.sample_id);
    const csv = formatDistanceMatrixCsv(ids, distMatrix, snpMatrix, sitesMatrix);
    triggerBrowserDownload(csv, `viromics_${selectedRefName || 'tree'}_distances.csv`, 'text/csv;charset=utf-8;');
  };

  const handleDownloadSvg = () => {
    if (!svgRef.current) return;
    const svgText = serializeSvgElement(svgRef.current);
    triggerBrowserDownload(svgText, `viromics_${selectedRefName || 'tree'}_${algorithm.toLowerCase()}.svg`, 'image/svg+xml;charset=utf-8;');
  };

  const handleDownloadCohortFasta = () => {
    const fasta = formatCohortMultiFasta(includedSamples, selectedRefName, { minBreadth });
    if (!fasta) return;
    triggerBrowserDownload(fasta, `viromics_${selectedRefName || 'cohort'}_consensus.fasta`, 'text/plain;charset=utf-8;');
  };

  const totalExcluded = excludedLowBreadth.length + excludedOverlap.length;

  return (
    <div className="p-5 glass-panel flex flex-col gap-4">
      <div className="flex flex-col xl:flex-row xl:items-center justify-between gap-3">
        <div>
          <h3 className="text-[16px] font-semibold text-[#111827]">
            Phylogenetic Tree ({algorithm === 'NJ' ? 'Midpoint-Rooted Neighbor-Joining' : 'Rooted UPGMA'})
          </h3>
          <p className="text-[12px] text-[#6B7280] mt-0.5">
            Coordinate-aligned, deletion-masked SNP divergence with zero-overlap guards and exact compared-site tracking
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-2.5 shrink-0">
          <div className="flex bg-[#F1F3F5] border border-[#E2E6EA] rounded-[4px] p-0.5 text-xs">
            <button
              onClick={() => setAlgorithm('NJ')}
              className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                algorithm === 'NJ' ? 'bg-[#0F5C6E] text-white font-semibold' : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              Neighbor-Joining (Midpoint)
            </button>
            <button
              onClick={() => setAlgorithm('UPGMA')}
              className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                algorithm === 'UPGMA' ? 'bg-[#0F5C6E] text-white font-semibold' : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              UPGMA
            </button>
          </div>

          <div className="flex items-center gap-1.5 bg-white border border-[#CBD2D9] rounded-[4px] px-2.5 py-1 text-xs text-[#4B5563]">
            <span>Min Breadth %:</span>
            <input
              type="number"
              min={0}
              max={100}
              step={5}
              value={minBreadth}
              onChange={(e) => setMinBreadth(Math.max(0, Math.min(100, Number(e.target.value) || 0)))}
              className="w-12 bg-[#F8F9FA] border border-[#CBD2D9] rounded-[3px] px-1 py-0.5 text-xs font-mono text-[#111827] tabular-nums"
            />
          </div>

          <div className="flex items-center gap-1.5 bg-white border border-[#CBD2D9] rounded-[4px] px-2.5 py-1 text-xs text-[#4B5563]">
            <span>Min Sites:</span>
            <input
              type="number"
              min={1}
              step={50}
              value={minComparedSites}
              onChange={(e) => setMinComparedSites(Math.max(1, Number(e.target.value) || 1))}
              className="w-16 bg-[#F8F9FA] border border-[#CBD2D9] rounded-[3px] px-1 py-0.5 text-xs font-mono text-[#111827] tabular-nums"
            />
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

      <div className="flex flex-wrap items-center justify-between gap-3 border-t border-b border-[#E2E6EA] py-2.5">
        <div className="flex items-center gap-3 text-xs text-[#4B5563] tabular-nums">
          <span>
            Included: <strong className="font-mono text-[#111827]">{includedSamples.length}</strong>
          </span>
          <span className="text-[#CBD2D9]">·</span>
          <span>
            Excluded: <strong className="font-mono text-[#111827]">{totalExcluded}</strong>
          </span>
        </div>

        {treeData && (
          <div className="flex flex-wrap items-center gap-2">
            <button
              onClick={handleDownloadNewick}
              className="flex items-center gap-1.5 px-2.5 py-1.5 bg-white border border-[#CBD2D9] text-[#111827] text-xs font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] transition-colors"
            >
              <Download className="w-3.5 h-3.5 text-[#0F5C6E]" />
              Newick (.nwk)
            </button>
            <button
              onClick={handleDownloadMatrixCsv}
              className="flex items-center gap-1.5 px-2.5 py-1.5 bg-white border border-[#CBD2D9] text-[#111827] text-xs font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] transition-colors"
            >
              <Download className="w-3.5 h-3.5 text-[#0F5C6E]" />
              Matrix (.csv)
            </button>
            <button
              onClick={handleDownloadSvg}
              className="flex items-center gap-1.5 px-2.5 py-1.5 bg-white border border-[#CBD2D9] text-[#111827] text-xs font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] transition-colors"
            >
              <Download className="w-3.5 h-3.5 text-[#0F5C6E]" />
              SVG
            </button>
            <button
              onClick={handleDownloadCohortFasta}
              className="flex items-center gap-1.5 px-3 py-1.5 bg-[#0F5C6E] hover:bg-[#0B4654] text-white text-xs font-semibold rounded-[4px] transition-colors"
            >
              <Download className="w-3.5 h-3.5" />
              Cohort FASTA
            </button>
          </div>
        )}
      </div>

      <div className="bg-[#F8F9FA] border border-[#E2E6EA] rounded-[4px] px-3.5 py-2 flex items-center gap-2.5 text-xs text-[#4B5563]">
        <Info className="w-4 h-4 text-[#0F5C6E] shrink-0" />
        <span>
          <strong className="text-[#111827]">Quality &amp; Overlap Safe-Guard (&ge;{minBreadth}% breadth, &ge;{minComparedSites.toLocaleString()} comparable sites):</strong>{' '}
          Excludes samples with insufficient mutual coverage or <code className="font-mono text-[#0B4654]">null</code> distance; negative NJ branch lengths are clamped to 0.
        </span>
      </div>

      {includedSamples.length < 2 ? (
        <div className="p-8 text-center text-[#6B7280] border border-[#E2E6EA] rounded-[6px] bg-[#F8F9FA] flex flex-col items-center gap-2">
          <AlertCircle className="w-6 h-6 text-[#8C5308]" />
          <p className="font-semibold text-[#111827]">Insufficient Comparable Samples</p>
          <p className="text-xs text-[#4B5563]">
            Tree reconstruction requires at least 2 mutually overlapping samples with &ge;{minBreadth}% coverage breadth for{' '}
            <span className="font-mono text-[#0F5C6E] font-semibold">{selectedRefName}</span> ({includedSamples.length} qualified).
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-5">
          <div className="lg:col-span-3 border border-[#E2E6EA] rounded-[6px] bg-white p-4 overflow-x-auto">
            {treeData && (
              <svg
                ref={svgRef}
                viewBox={`0 0 ${layout.svgWidth} ${layout.svgHeight}`}
                width="100%"
                className="min-w-[620px] select-none"
                style={{ fontFamily: 'Arial, Helvetica, sans-serif' }}
              >
                {[0, 0.25, 0.5, 0.75, 1].map((pct, i) => {
                  const x = layout.leftMargin + pct * layout.treeWidth;
                  const scaleVal =
                    algorithm === 'UPGMA'
                      ? pct * layout.maxDepth * 2 * 100
                      : pct * layout.maxDepth * 100;
                  return (
                    <g key={i}>
                      <line
                        x1={x}
                        y1={16}
                        x2={x}
                        y2={layout.svgHeight - 28}
                        stroke="#E5E7EB"
                        strokeWidth="1"
                        strokeDasharray="3 3"
                      />
                      <text
                        x={x}
                        y={layout.svgHeight - 10}
                        fill="#6B7280"
                        fontSize="10"
                        textAnchor="middle"
                      >
                        {scaleVal.toFixed(3)}%
                      </text>
                    </g>
                  );
                })}

                {treeEdges.map((edge) => (
                  <path
                    key={edge.key}
                    d={edge.d}
                    fill="none"
                    stroke="#2D3742"
                    strokeWidth="1.5"
                    strokeLinecap="round"
                  />
                ))}

                {Object.values(layout.coords).map((c) => {
                  const isLeaf = c.node.isLeaf;
                  const isHovered = hoveredNode && hoveredNode.id === c.node.id;
                  return (
                    <g
                      key={c.node.id}
                      className="cursor-pointer"
                      onMouseEnter={() => setHoveredNode(c.node)}
                      onMouseLeave={() => setHoveredNode(null)}
                    >
                      <circle
                        cx={c.x}
                        cy={c.y}
                        r={isLeaf ? 4 : 3}
                        fill={isHovered ? '#B43B2E' : isLeaf ? '#0F5C6E' : '#6B7280'}
                      />

                      {isLeaf ? (
                        <text
                          x={c.x + 10}
                          y={c.y + 4}
                          fill={isHovered ? '#B43B2E' : '#111827'}
                          fontSize="11.5"
                          fontWeight={isHovered ? '700' : '600'}
                          textAnchor="start"
                        >
                          {c.node.name}
                        </text>
                      ) : (
                        isHovered && (
                          <g>
                            <rect
                              x={c.x - 52}
                              y={c.y - 24}
                              width="104"
                              height="17"
                              rx="4"
                              fill="#FFFFFF"
                              stroke="#CBD2D9"
                              strokeWidth="1"
                            />
                            <text
                              x={c.x}
                              y={c.y - 12.5}
                              fill="#111827"
                              fontSize="9"
                              fontWeight="600"
                              textAnchor="middle"
                            >
                              {`Len: ${((c.node.branchLength || 0) * 100).toFixed(3)}%`}
                            </text>
                          </g>
                        )
                      )}
                    </g>
                  );
                })}
              </svg>
            )}
          </div>

          <div className="flex flex-col gap-4">
            <div className="p-4 bg-[#F8F9FA] border border-[#E2E6EA] rounded-[6px]">
              <div className="flex items-center gap-2 text-[#0F5C6E] mb-1.5">
                <GitBranch className="w-4 h-4" />
                <h4 className="text-[11px] font-semibold tracking-[0.06em] uppercase">Nearest Relatives</h4>
              </div>
              <p className="text-[11px] text-[#6B7280] leading-relaxed">
                Pairwise SNP divergence, raw SNPs, and mutually covered comparable sites:
              </p>

              <div className="mt-3 flex flex-col gap-2.5 max-h-72 overflow-y-auto pr-1">
                {includedSamples.map((s, idx) => {
                  let minVal = Infinity;
                  let matchIdx = -1;
                  (distMatrix[idx] || []).forEach((val, j) => {
                    if (idx !== j && val !== null && val < minVal && includedSamples[j]) {
                      minVal = val;
                      matchIdx = j;
                    }
                  });
                  const snpVal = matchIdx >= 0 ? snpMatrix?.[idx]?.[matchIdx] : null;
                  const sitesVal = matchIdx >= 0 ? sitesMatrix?.[idx]?.[matchIdx] : null;

                  return (
                    <div key={s.sample_id} className="flex flex-col border-b border-[#E2E6EA] pb-2 last:border-none">
                      <span className="text-[12px] font-semibold text-[#111827] truncate">{s.sample_id}</span>
                      {matchIdx >= 0 ? (
                        <div className="text-[11px] text-[#4B5563] mt-0.5 space-y-0.5 tabular-nums">
                          <div>
                            Nearest: <strong className="text-[#0F5C6E] font-semibold">{includedSamples[matchIdx].sample_id}</strong>
                          </div>
                          <div className="font-mono text-[#4B5563]">
                            {(minVal * 100).toFixed(3)}% div
                            {snpVal !== null && snpVal !== undefined ? ` · ${snpVal} SNP${snpVal === 1 ? '' : 's'}` : ''}
                            {sitesVal !== null && sitesVal !== undefined ? ` · ${sitesVal.toLocaleString()} bp` : ''}
                          </div>
                        </div>
                      ) : (
                        <span className="text-[11px] text-[#6B7280] mt-0.5">No comparable relative</span>
                      )}
                    </div>
                  );
                })}
              </div>
            </div>
          </div>
        </div>
      )}

      {(excludedLowBreadth.length > 0 || excludedOverlap.length > 0) && (
        <div className="border border-[#E2E6EA] bg-[#F8F9FA] rounded-[4px] px-4 py-3 flex flex-col gap-1.5">
          <span className="text-xs font-semibold text-[#A62C22]">
            Excluded Samples (Low Breadth &lt;{minBreadth}% at &ge;{minCoverageThreshold}x or Insufficient Overlap)
          </span>
          <div className="flex flex-wrap gap-2">
            {excludedLowBreadth.map(({ sample, breadth }) => (
              <span key={`lb-${sample.sample_id}`} className="text-[11px] bg-white border border-[#E2E6EA] text-[#111827] font-mono px-2.5 py-1 rounded-[4px] tabular-nums">
                {sample.sample_id}: <strong className="text-[#A62C22]">{breadth.toFixed(1)}% cov</strong>
              </span>
            ))}
            {excludedOverlap.map(({ sample, breadth }) => (
              <span key={`ov-${sample.sample_id}`} className="text-[11px] bg-white border border-[#EBD29B] text-[#111827] font-mono px-2.5 py-1 rounded-[4px] tabular-nums">
                {sample.sample_id}: <strong className="text-[#8C5308]">{breadth.toFixed(1)}% cov (insufficient site overlap)</strong>
              </span>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
