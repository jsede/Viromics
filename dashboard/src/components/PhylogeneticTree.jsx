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
      <div className="p-8 text-center text-slate-500 glass-panel rounded-2xl flex flex-col items-center gap-3">
        <AlertCircle className="w-8 h-8 text-amber-500" />
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

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      <div className="flex flex-col xl:flex-row xl:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold text-slate-900">
            Phylogenetic Tree ({algorithm === 'NJ' ? 'Midpoint-Rooted Neighbor-Joining' : 'Rooted UPGMA'})
          </h3>
          <p className="text-xs text-slate-500">
            Built from coordinate-aligned, deletion-masked SNP divergence with zero-overlap guards and exact compared-site tracking
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-2.5">
          <div className="flex bg-slate-100 border border-slate-200 rounded-lg p-0.5 text-xs">
            <button
              onClick={() => setAlgorithm('NJ')}
              className={`px-3 py-1 rounded-md font-medium transition-colors ${
                algorithm === 'NJ' ? 'bg-red-600 text-white font-semibold' : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              Neighbor-Joining (Midpoint)
            </button>
            <button
              onClick={() => setAlgorithm('UPGMA')}
              className={`px-3 py-1 rounded-md font-medium transition-colors ${
                algorithm === 'UPGMA' ? 'bg-red-600 text-white font-semibold' : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              UPGMA
            </button>
          </div>

          <div className="flex items-center gap-1.5 bg-slate-50 border border-slate-200 rounded-xl px-2.5 py-1 text-xs text-slate-700">
            <span>Min Breadth %:</span>
            <input
              type="number"
              min={0}
              max={100}
              step={5}
              value={minBreadth}
              onChange={(e) => setMinBreadth(Math.max(0, Math.min(100, Number(e.target.value) || 0)))}
              className="w-12 bg-white border border-slate-300 rounded px-1 py-0.5 text-xs font-mono text-slate-900"
            />
          </div>

          <div className="flex items-center gap-1.5 bg-slate-50 border border-slate-200 rounded-xl px-2.5 py-1 text-xs text-slate-700">
            <span>Min Sites:</span>
            <input
              type="number"
              min={1}
              step={50}
              value={minComparedSites}
              onChange={(e) => setMinComparedSites(Math.max(1, Number(e.target.value) || 1))}
              className="w-16 bg-white border border-slate-300 rounded px-1 py-0.5 text-xs font-mono text-slate-900"
            />
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

      <div className="bg-indigo-50/70 border border-indigo-200 rounded-xl p-3.5 flex flex-col sm:flex-row sm:items-center justify-between gap-3">
        <div className="flex items-start gap-3">
          <Info className="w-5 h-5 text-indigo-600 shrink-0 mt-0.5" />
          <div>
            <span className="text-xs font-semibold text-indigo-900 block">
              Quality &amp; Overlap Safe-Guard (&ge;{minBreadth}% breadth, &ge;{minComparedSites.toLocaleString()} comparable sites)
            </span>
            <span className="text-[11px] text-slate-600 block mt-0.5 leading-relaxed">
              Excludes samples below {minBreadth}% breadth or with <code>null</code> / insufficient mutual coverage to prevent topological distortion. Negative NJ branch lengths are clamped to 0.
            </span>
          </div>
        </div>

        {treeData && (
          <div className="flex flex-wrap items-center gap-2 shrink-0">
            <button
              onClick={handleDownloadNewick}
              className="flex items-center gap-1 px-2.5 py-1.5 bg-white border border-slate-300 text-slate-700 text-[11px] font-semibold rounded-lg hover:bg-slate-50"
            >
              <Download className="w-3.5 h-3.5 text-red-600" />
              Newick (.nwk)
            </button>
            <button
              onClick={handleDownloadMatrixCsv}
              className="flex items-center gap-1 px-2.5 py-1.5 bg-white border border-slate-300 text-slate-700 text-[11px] font-semibold rounded-lg hover:bg-slate-50"
            >
              <Download className="w-3.5 h-3.5 text-red-600" />
              Matrix (.csv)
            </button>
            <button
              onClick={handleDownloadSvg}
              className="flex items-center gap-1 px-2.5 py-1.5 bg-white border border-slate-300 text-slate-700 text-[11px] font-semibold rounded-lg hover:bg-slate-50"
            >
              <Download className="w-3.5 h-3.5 text-red-600" />
              SVG
            </button>
            <button
              onClick={handleDownloadCohortFasta}
              className="flex items-center gap-1 px-2.5 py-1.5 bg-red-600 text-white text-[11px] font-semibold rounded-lg hover:bg-red-700"
            >
              <Download className="w-3.5 h-3.5" />
              Cohort FASTA
            </button>
          </div>
        )}
      </div>

      {includedSamples.length < 2 ? (
        <div className="p-8 text-center text-slate-500 border border-slate-200 rounded-xl bg-slate-50 flex flex-col items-center gap-2">
          <AlertCircle className="w-8 h-8 text-amber-500" />
          <p className="font-semibold text-slate-900">Insufficient Comparable Samples</p>
          <p className="text-xs text-slate-600">
            Tree reconstruction requires at least 2 mutually overlapping samples with &ge;{minBreadth}% coverage breadth for{' '}
            <span className="font-mono text-red-600 font-semibold">{selectedRefName}</span> ({includedSamples.length} qualified).
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          <div className="lg:col-span-3 border border-slate-200 rounded-xl bg-white p-4 overflow-x-auto">
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
                        stroke="#E2E8F0"
                        strokeWidth="1"
                        strokeDasharray="3 3"
                      />
                      <text
                        x={x}
                        y={layout.svgHeight - 10}
                        fill="#64748B"
                        fontSize="9.5"
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
                    stroke="#475569"
                    strokeWidth="2"
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
                        r={isLeaf ? 4.5 : 3.5}
                        fill={isHovered ? '#DC2626' : isLeaf ? '#0F766E' : '#4F46E5'}
                      />

                      {isLeaf ? (
                        <text
                          x={c.x + 10}
                          y={c.y + 4}
                          fill={isHovered ? '#DC2626' : '#0F172A'}
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
                              height="16"
                              rx="4"
                              fill="#0F172A"
                            />
                            <text
                              x={c.x}
                              y={c.y - 13}
                              fill="#FFFFFF"
                              fontSize="8.5"
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
            <div className="p-4 bg-slate-50 border border-slate-200 rounded-xl">
              <div className="flex items-center gap-2 text-red-600 mb-2">
                <GitBranch className="w-4 h-4" />
                <h4 className="text-xs font-bold tracking-wider uppercase">Nearest Relatives</h4>
              </div>
              <p className="text-[11px] text-slate-600 leading-relaxed">
                Pairwise SNP divergence, raw SNP differences, and mutually covered comparable sites:
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
                    <div key={s.sample_id} className="flex flex-col border-b border-slate-200 pb-2 last:border-none">
                      <span className="text-[11px] font-semibold text-slate-900 truncate">{s.sample_id}</span>
                      {matchIdx >= 0 ? (
                        <div className="text-[10px] text-slate-600 mt-0.5 space-y-0.5">
                          <div>
                            Nearest: <strong className="text-red-600">{includedSamples[matchIdx].sample_id}</strong>
                          </div>
                          <div className="font-mono text-slate-700">
                            {(minVal * 100).toFixed(3)}% div
                            {snpVal !== null && snpVal !== undefined ? ` · ${snpVal} SNP${snpVal === 1 ? '' : 's'}` : ''}
                            {sitesVal !== null && sitesVal !== undefined ? ` · ${sitesVal.toLocaleString()} bp` : ''}
                          </div>
                        </div>
                      ) : (
                        <span className="text-[10px] text-slate-400 mt-0.5">No comparable relative</span>
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
        <div className="border border-slate-200 bg-slate-50 rounded-xl p-4 flex flex-col gap-2">
          <span className="text-xs font-semibold text-red-700 tracking-wider block">
            Excluded Samples (Low Breadth &lt;{minBreadth}% at &ge;{minCoverageThreshold}x or Insufficient Overlap)
          </span>
          <div className="flex flex-wrap gap-2 mt-1">
            {excludedLowBreadth.map(({ sample, breadth }) => (
              <span key={`lb-${sample.sample_id}`} className="text-[11px] bg-white border border-slate-200 text-slate-700 font-mono px-2.5 py-1 rounded-lg">
                {sample.sample_id}: <strong className="text-red-600">{breadth.toFixed(1)}% cov</strong>
              </span>
            ))}
            {excludedOverlap.map(({ sample, breadth }) => (
              <span key={`ov-${sample.sample_id}`} className="text-[11px] bg-white border border-amber-200 text-slate-700 font-mono px-2.5 py-1 rounded-lg">
                {sample.sample_id}: <strong className="text-amber-600">{breadth.toFixed(1)}% cov (insufficient site overlap)</strong>
              </span>
            ))}
          </div>
        </div>
      )}
    </div>
  );
}
