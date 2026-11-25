import React, { useState, useEffect } from 'react';
import { GitBranch, AlertCircle, Info, Download } from 'lucide-react';

export default function PhylogeneticTree({ 
  samples, 
  selectedRefName, 
  setSelectedRefName, 
  refs = [], 
  distanceMatrices = {}, 
  distanceMatrixOrder = [],
  minCoverageThreshold = 10
}) {
  const [treeData, setTreeData] = useState(null);
  const [distanceMatrix, setDistanceMatrix] = useState([]);
  const [hoveredNode, setHoveredNode] = useState(null);

  const getBreadth = (sample) => {
    const cov = (selectedRefName && sample.coverage && sample.coverage[selectedRefName]) ? sample.coverage[selectedRefName] : (sample.coverage || {});
    return cov.coverage_breadth_min_cov !== undefined ? cov.coverage_breadth_min_cov : (cov.coverage_breadth_10x !== undefined ? cov.coverage_breadth_10x : 0);
  };

  // Filter samples into passing (>= 80% coverage breadth) and excluded (< 80%) for the active strain
  const passingSamples = samples.filter(sample => getBreadth(sample) >= 80.0);
  const excludedSamples = samples.filter(sample => getBreadth(sample) < 80.0);

  useEffect(() => {
    if (!samples || passingSamples.length < 2) {
      setTreeData(null);
      setDistanceMatrix([]);
      return;
    }
    
    // 1. Load and slice the pre-calculated distance matrix for passing samples
    const numSamples = passingSamples.length;
    const matrix = Array(numSamples).fill(0).map(() => Array(numSamples).fill(0));
    
    // Get the pre-calculated matrix for the selected reference
    const precalcMatrix = (selectedRefName && distanceMatrices && distanceMatrices[selectedRefName]) ? distanceMatrices[selectedRefName] : null;
    
    if (precalcMatrix) {
      // Find the index of each passing sample in the distance matrix order (if available) to ensure stable alignment
      const originalIndices = passingSamples.map(ps => {
        if (distanceMatrixOrder && Array.isArray(distanceMatrixOrder) && distanceMatrixOrder.length > 0) {
          return distanceMatrixOrder.findIndex(name => name === ps.sample_id);
        }
        return samples.findIndex(s => s.sample_id === ps.sample_id);
      });
      
      for (let i = 0; i < numSamples; i++) {
        const idxA = originalIndices[i];
        for (let j = 0; j < numSamples; j++) {
          const idxB = originalIndices[j];
          if (idxA !== -1 && idxB !== -1 && precalcMatrix[idxA] && precalcMatrix[idxA][idxB] !== undefined) {
            matrix[i][j] = precalcMatrix[idxA][idxB];
          } else {
            matrix[i][j] = 0.05; // Fallback
          }
        }
      }
    } else {
      // Fallback distance calculation if no precalculated matrix exists (e.g. legacy/mock reports)
      const getSeq = (sample) => {
        if (!sample || !sample.consensus) return null;
        const c = (selectedRefName && sample.consensus[selectedRefName]) ? sample.consensus[selectedRefName] : sample.consensus;
        if (!c) return null;
        if (typeof c === 'string') return c;
        if (typeof c.sequence === 'string') return c.sequence;
        return null;
      };

      for (let i = 0; i < numSamples; i++) {
        for (let j = i + 1; j < numSamples; j++) {
          const seqA = getSeq(passingSamples[i]);
          const seqB = getSeq(passingSamples[j]);
          
          let dist = 0.05;
          if (seqA && seqB) {
            let diffs = 0;
            let compared = 0;
            const len = Math.min(seqA.length, seqB.length);
            for (let k = 0; k < len; k++) {
              const charA = seqA[k];
              const charB = seqB[k];
              if (charA !== 'N' && charB !== 'N') {
                compared++;
                if (charA !== charB) {
                  diffs++;
                }
              }
            }
            dist = compared > 0 ? (diffs / compared) : 0;
          }
          matrix[i][j] = dist;
          matrix[j][i] = dist;
        }
      }
    }
    
    setDistanceMatrix(matrix);

    // 2. Perform UPGMA Hierarchical Clustering to build rooted tree
    // Initialize clusters (nodes)
    let clusters = passingSamples.map((s, idx) => ({
      id: `leaf_${idx}`,
      name: s.sample_id,
      isLeaf: true,
      height: 0,
      indices: [idx],
      children: []
    }));

    // Keep active clusters and distances
    let activeClusters = [...clusters];

    while (activeClusters.length > 1) {
      // Find pair of clusters with minimum distance
      let minVal = Infinity;
      let u = -1, v = -1;
      
      for (let i = 0; i < activeClusters.length; i++) {
        for (let j = i + 1; j < activeClusters.length; j++) {
          // Calculate distance between cluster i and cluster j (average distance)
          let distSum = 0;
          let count = 0;
          
          activeClusters[i].indices.forEach(idxA => {
            activeClusters[j].indices.forEach(idxB => {
              distSum += matrix[idxA][idxB];
              count++;
            });
          });
          
          const avgDist = count > 0 ? (distSum / count) : 0.05;
          if (avgDist < minVal) {
            minVal = avgDist;
            u = i;
            v = j;
          }
        }
      }

      if (u === -1 || v === -1) break;

      // Merge clusters u and v
      const nodeA = activeClusters[u];
      const nodeB = activeClusters[v];
      const newHeight = minVal / 2;

      const mergedNode = {
        id: `node_${Date.now()}_${Math.random()}`,
        name: `Dist: ${(minVal * 100).toFixed(3)}%`,
        isLeaf: false,
        height: newHeight,
        indices: [...nodeA.indices, ...nodeB.indices],
        children: [nodeA, nodeB]
      };

      // Remove u and v (v > u, so remove v first to preserve index)
      activeClusters.splice(v, 1);
      activeClusters.splice(u, 1);
      
      // Add merged cluster
      activeClusters.push(mergedNode);
    }

    if (activeClusters.length > 0) {
      setTreeData(activeClusters[0]);
    }
  }, [samples, selectedRefName, distanceMatrices, distanceMatrixOrder]);

  if (!samples || samples.length < 2) {
    return (
      <div className="p-8 text-center text-gray-400 glass-panel rounded-2xl flex flex-col items-center gap-3">
        <AlertCircle className="w-8 h-8 text-brand-amber animate-pulse" />
        <p>Phylogenetic clustering requires at least 2 samples in the run.</p>
      </div>
    );
  }

  // Linear layout function for the tree visualization
  const getCoordinates = (node, totalLeaves, maxTreeHeight) => {
    let leavesCount = 0;
    const coords = {};

    const assignCoords = (n, xStart, depth) => {
      const x = 50 + (1 - (maxTreeHeight > 0 ? (n.height / maxTreeHeight) : 0)) * 600; // tree scaled to width=600, starting x=50
      
      if (n.isLeaf) {
        const y = 40 + leavesCount * 50; // leaves spaced vertically by 50px
        leavesCount++;
        coords[n.id] = { x, y, node: n };
        return y;
      } else {
        const yCoords = n.children.map(c => assignCoords(c, x, depth + 1));
        const y = yCoords.reduce((a, b) => a + b, 0) / yCoords.length; // center parent between children
        coords[n.id] = { x, y, node: n };
        return y;
      }
    };

    if (node) {
      assignCoords(node, 50, 0);
    }
    return coords;
  };

  const maxTreeHeight = treeData ? treeData.height : 0.05;
  const numLeaves = passingSamples.length;
  const nodeCoordinates = treeData ? getCoordinates(treeData, numLeaves, maxTreeHeight) : {};
  const svgHeight = Math.max(250, numLeaves * 60);

  // SVG lines generation helper
  const drawTreeLines = (node, coords, lines = []) => {
    if (node && !node.isLeaf && node.children.length === 2) {
      const parent = coords[node.id];
      const childA = coords[node.children[0].id];
      const childB = coords[node.children[1].id];

      if (parent && childA && childB) {
        // Draw UPGMA square tree lines (horizontal to child X, vertical connecting children)
        lines.push(
          <path
            key={`line_${node.id}`}
            d={`M ${parent.x} ${parent.y} L ${parent.x} ${childA.y} L ${childA.x} ${childA.y} M ${parent.x} ${parent.y} L ${parent.x} ${childB.y} L ${childB.x} ${childB.y}`}
            fill="none"
            stroke="rgba(99, 102, 241, 0.6)"
            strokeWidth="2.5"
            strokeLinecap="round"
          />
        );
      }
      drawTreeLines(node.children[0], coords, lines);
      drawTreeLines(node.children[1], coords, lines);
    }
    return lines;
  };

  const treePaths = treeData ? drawTreeLines(treeData, nodeCoordinates) : [];

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      {/* Header Row with dynamic dropdown */}
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold font-outfit text-white">Consensus Distance Tree</h3>
          <p className="text-xs text-gray-400">
            Rooted UPGMA dendrogram calculated in the browser based on consensus sequence pairwise divergence (Hamming Distance).
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

      {/* Informative Warning Note on Coverage cutoff */}
      <div className="bg-indigo-950/20 border border-brand-indigo/20 rounded-xl p-3.5 flex items-start gap-3">
        <Info className="w-5 h-5 text-brand-indigo shrink-0 mt-0.5" />
        <div className="flex-1">
          <span className="text-xs font-semibold text-brand-indigo block font-outfit">
            Epidemiological Quality Safe-Guard (≥80% Coverage Required)
          </span>
          <span className="text-[11px] text-gray-400 block mt-1 leading-relaxed">
            To prevent false-positive clustering and branch distortion caused by sequencing dropouts (long stretches of uncalled NNNs), this distance tree strictly filters out any sample with <strong>less than 80% coverage breadth</strong> for the active strain. The Mutational Sharing Matrix continues to display all samples above 50% coverage for broader mutation tracking.
          </span>
        </div>
      </div>

      {/* Main Tree & Matrix Content */}
      {passingSamples.length < 2 ? (
        <div className="p-8 text-center text-gray-400 border border-dark-border rounded-xl bg-dark-bg/20 flex flex-col items-center gap-3">
          <AlertCircle className="w-8 h-8 text-brand-amber animate-pulse" />
          <p className="font-semibold text-white">Insufficient Samples with Coverage</p>
          <p className="text-xs text-gray-400">
            Phylogenetic clustering requires at least 2 samples with ≥ 80% coverage breadth for reference: <span className="font-mono text-brand-teal">{selectedRefName}</span>.
            Currently only {passingSamples.length} sample(s) passed the threshold.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-4 gap-6">
          {/* SVG Drawing Canvas */}
          <div className="lg:col-span-3 border border-dark-border rounded-xl bg-dark-bg/40 p-4 overflow-x-auto">
            {treeData && (
              <svg viewBox={`0 0 800 ${svgHeight}`} width="100%" className="tree-svg min-w-[600px] select-none">
                {/* Scale grid */}
                {[0, 0.25, 0.5, 0.75, 1].map((pct, i) => {
                  const dist = pct * maxTreeHeight;
                  const x = 50 + (1 - pct) * 600;
                  return (
                    <g key={i}>
                      <line x1={x} y1={20} x2={x} y2={svgHeight - 20} stroke="rgba(36,47,77,0.2)" strokeDasharray="3 3" />
                      <text x={x} y={svgHeight - 5} fill="#9ca3af" fontSize="8" textAnchor="middle" fontFamily="Fira Code">
                        {(dist * 2 * 100).toFixed(2)}% div
                      </text>
                    </g>
                  );
                })}

                {/* Tree lines */}
                {treePaths}

                {/* Tree nodes (labels + circles) */}
                {Object.values(nodeCoordinates).map((c, i) => {
                  const isLeaf = c.node.isLeaf;
                  const isHovered = hoveredNode && hoveredNode.id === c.node.id;

                  return (
                    <g 
                      key={i} 
                      className="cursor-pointer"
                      onMouseEnter={() => setHoveredNode(c.node)}
                      onMouseLeave={() => setHoveredNode(null)}
                    >
                      <circle 
                        cx={c.x} cy={c.y} 
                        r={isLeaf ? 5 : 4} 
                        className={`transition-all duration-150 ${
                          isLeaf 
                            ? isHovered ? 'fill-brand-rose' : 'fill-brand-teal' 
                            : isHovered ? 'fill-brand-rose' : 'fill-brand-indigo'
                        }`}
                      />
                      
                      {isLeaf ? (
                        <text 
                          x={c.x + 12} y={c.y + 4}
                          fill={isHovered ? '#fff' : '#e5e7eb'}
                          fontSize="11" fontWeight={isHovered ? '700' : '500'}
                          fontFamily="Outfit"
                          className="transition-colors"
                        >
                          {c.node.name}
                        </text>
                      ) : (
                        isHovered && (
                          <g>
                            <rect 
                              x={c.x - 30} y={c.y - 22} width="60" height="15" rx="3" 
                              fill="#161d30" stroke="#242f4d" strokeWidth="1"
                            />
                            <text 
                              x={c.x} y={c.y - 12}
                              fill="#8b5cf6" fontSize="8" fontWeight="bold" textAnchor="middle"
                              fontFamily="Fira Code"
                            >
                              {c.node.name}
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

          {/* Distance Matrix Summary */}
          <div className="flex flex-col gap-4">
            <div className="p-4 bg-dark-bg/60 border border-dark-border rounded-xl">
              <div className="flex items-center gap-2 text-brand-teal mb-2">
                <GitBranch className="w-4 h-4" />
                <h4 className="text-xs font-bold font-outfit uppercase tracking-wider">Distance Matrix</h4>
              </div>
              <p className="text-[11px] text-gray-400 leading-relaxed">
                Pairwise genetic divergence matches the number of nucleotide mutations compared over the non-N genome sequence length.
              </p>
              
              {/* Quick list of sample relationships */}
              <div className="mt-4 flex flex-col gap-2.5 max-h-48 overflow-y-auto">
                {passingSamples.map((s, idx) => {
                  // Find closest match
                  let minVal = Infinity;
                  let matchId = null;
                  distanceMatrix[idx]?.forEach((val, j) => {
                    if (idx !== j && val < minVal && passingSamples[j]) {
                      minVal = val;
                      matchId = passingSamples[j].sample_id;
                    }
                  });

                  return (
                    <div key={idx} className="flex flex-col border-b border-dark-border/40 pb-2">
                      <span className="text-[11px] font-semibold text-gray-200 truncate">{s.sample_id}</span>
                      {matchId ? (
                        <span className="text-[10px] text-gray-400 mt-0.5">
                          Nearest: <strong className="text-brand-teal truncate inline-block max-w-[120px] align-bottom">{matchId}</strong>{' '}
                          (<span className="font-mono text-brand-teal">{(minVal*100).toFixed(3)}%</span>)
                        </span>
                      ) : (
                        <span className="text-[10px] text-gray-500 mt-0.5">No relative found</span>
                      )}
                    </div>
                  );
                })}
              </div>
            </div>
          </div>
        </div>
      )}

      {/* Excluded Samples Box */}
      {excludedSamples.length > 0 && (
        <div className="border border-dark-border/60 bg-dark-bg/20 rounded-xl p-4 flex flex-col gap-2">
          <span className="text-xs font-semibold text-brand-rose uppercase tracking-wider font-outfit block">
            Excluded Low-Coverage Samples (&lt;80% Breadth at ≥{minCoverageThreshold}x)
          </span>
          <p className="text-[11px] text-gray-400">
            The following samples had insufficient coverage breadth for this strain and were excluded from phylogenetic distance analysis to prevent artificial clustering on uncalled regions (NNNs):
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
