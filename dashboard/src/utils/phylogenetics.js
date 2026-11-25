import { getEvalBreadth, getContigConsensus } from './reportSchema.js';

/**
 * Filter candidate samples by breadth and pairwise overlap / minComparedSites.
 * Guarantees a complete numeric submatrix with zero nulls for tree reconstruction.
 */
export function filterTreeSamples(
  samples,
  refName,
  report = {},
  { minBreadth = 80, minComparedSites = 1 } = {}
) {
  const safeSamples = Array.isArray(samples) ? samples : [];
  const excludedLowBreadth = [];
  let candidates = [];

  safeSamples.forEach((s) => {
    const br = getEvalBreadth(s, refName);
    if (br >= minBreadth) {
      candidates.push(s);
    } else {
      excludedLowBreadth.push({ sample: s, breadth: br, reason: 'low_breadth' });
    }
  });

  if (candidates.length < 2) {
    return {
      includedSamples: candidates,
      excludedLowBreadth,
      excludedOverlap: [],
      distMatrix: [],
      snpMatrix: [],
      sitesMatrix: [],
    };
  }

  const order = Array.isArray(report?.distance_matrix_order) ? report.distance_matrix_order : [];
  const preDist = refName ? report?.distance_matrices?.[refName] : null;
  const preSnp = refName ? report?.snp_distance_matrices?.[refName] : null;
  const preSites = refName ? report?.compared_sites_matrices?.[refName] : null;

  const getIndex = (sampleId) => {
    if (order.length > 0) return order.indexOf(sampleId);
    return safeSamples.findIndex((s) => s.sample_id === sampleId);
  };

  const rawPairData = (cands) => {
    const n = cands.length;
    const dMat = Array.from({ length: n }, () => Array(n).fill(null));
    const sMat = Array.from({ length: n }, () => Array(n).fill(null));
    const cMat = Array.from({ length: n }, () => Array(n).fill(null));

    if (preDist) {
      const idxs = cands.map((c) => getIndex(c.sample_id));
      for (let i = 0; i < n; i++) {
        const ia = idxs[i];
        dMat[i][i] = 0;
        sMat[i][i] = preSnp?.[ia]?.[ia] ?? 0;
        cMat[i][i] = preSites?.[ia]?.[ia] ?? null;
        for (let j = i + 1; j < n; j++) {
          const ib = idxs[j];
          const dv = ia !== -1 && ib !== -1 ? preDist?.[ia]?.[ib] ?? null : null;
          const sv = ia !== -1 && ib !== -1 && preSnp ? preSnp?.[ia]?.[ib] ?? null : null;
          const cv = ia !== -1 && ib !== -1 && preSites ? preSites?.[ia]?.[ib] ?? null : null;
          dMat[i][j] = dMat[j][i] = dv;
          sMat[i][j] = sMat[j][i] = sv;
          cMat[i][j] = cMat[j][i] = cv;
        }
      }
      return { dMat, sMat, cMat };
    }

    // Legacy fallback if report has no distance_matrices: consensus pairwise comparison (null on zero overlap)
    for (let i = 0; i < n; i++) {
      const seqA = getContigConsensus(cands[i], refName).sequence || '';
      const nonN = [...seqA].filter((ch) => ch !== 'N').length;
      dMat[i][i] = 0;
      sMat[i][i] = 0;
      cMat[i][i] = nonN;
      for (let j = i + 1; j < n; j++) {
        const seqB = getContigConsensus(cands[j], refName).sequence || '';
        const len = Math.min(seqA.length, seqB.length);
        let diffs = 0;
        let compared = 0;
        for (let k = 0; k < len; k++) {
          const a = seqA[k];
          const b = seqB[k];
          if (a !== 'N' && b !== 'N') {
            compared++;
            if (a.toUpperCase() !== b.toUpperCase()) diffs++;
          }
        }
        const dist = compared > 0 ? diffs / compared : null;
        dMat[i][j] = dMat[j][i] = dist;
        sMat[i][j] = sMat[j][i] = compared > 0 ? diffs : null;
        cMat[i][j] = cMat[j][i] = compared;
      }
    }
    return { dMat, sMat, cMat };
  };

  const excludedOverlap = [];
  let currentData = rawPairData(candidates);

  while (candidates.length >= 2) {
    const n = candidates.length;
    const badCounts = Array(n).fill(0);
    let hasBad = false;

    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const dv = currentData.dMat[i][j];
        const cv = currentData.cMat[i][j];
        const invalidDist = dv === null || dv === undefined || Number.isNaN(dv);
        const lowSites = cv !== null && cv !== undefined && cv < minComparedSites;
        if (invalidDist || lowSites) {
          badCounts[i]++;
          badCounts[j]++;
          hasBad = true;
        }
      }
    }

    if (!hasBad) break;

    let dropIdx = 0;
    for (let i = 1; i < n; i++) {
      if (badCounts[i] > badCounts[dropIdx]) {
        dropIdx = i;
      } else if (badCounts[i] === badCounts[dropIdx]) {
        const brI = getEvalBreadth(candidates[i], refName);
        const brDrop = getEvalBreadth(candidates[dropIdx], refName);
        if (brI < brDrop) dropIdx = i;
      }
    }

    const dropped = candidates[dropIdx];
    excludedOverlap.push({
      sample: dropped,
      breadth: getEvalBreadth(dropped, refName),
      reason: 'zero_overlap_or_insufficient_sites',
    });
    candidates = candidates.filter((_, idx) => idx !== dropIdx);
    currentData = candidates.length >= 2 ? rawPairData(candidates) : { dMat: [], sMat: [], cMat: [] };
  }

  return {
    includedSamples: candidates,
    excludedLowBreadth,
    excludedOverlap,
    distMatrix: currentData.dMat,
    snpMatrix: currentData.sMat,
    sitesMatrix: currentData.cMat,
  };
}

/**
 * Deterministic UPGMA hierarchical clustering returning tree with depthFromRoot & branchLength.
 */
export function buildUPGMATree(sampleIds, matrix) {
  if (!Array.isArray(sampleIds) || sampleIds.length === 0) return null;
  if (sampleIds.length === 1) {
    return {
      id: 'leaf_0',
      name: sampleIds[0],
      isLeaf: true,
      height: 0,
      branchLength: 0,
      depthFromRoot: 0,
      pairwiseDist: 0,
      indices: [0],
      children: [],
    };
  }

  let nextId = 1;
  let active = sampleIds.map((id, idx) => ({
    id: `leaf_${idx}`,
    name: id,
    isLeaf: true,
    height: 0,
    branchLength: 0,
    depthFromRoot: 0,
    pairwiseDist: 0,
    indices: [idx],
    children: [],
  }));

  while (active.length > 1) {
    let minVal = Infinity;
    let u = -1;
    let v = -1;

    for (let i = 0; i < active.length; i++) {
      for (let j = i + 1; j < active.length; j++) {
        let sum = 0;
        let cnt = 0;
        active[i].indices.forEach((ia) => {
          active[j].indices.forEach((ib) => {
            sum += matrix[ia][ib];
            cnt++;
          });
        });
        const avg = cnt > 0 ? sum / cnt : 0;
        if (avg < minVal) {
          minVal = avg;
          u = i;
          v = j;
        }
      }
    }

    if (u === -1 || v === -1) break;
    const nodeA = active[u];
    const nodeB = active[v];
    const newHeight = Math.max(0, minVal / 2);

    nodeA.branchLength = Math.max(0, newHeight - nodeA.height);
    nodeB.branchLength = Math.max(0, newHeight - nodeB.height);

    const merged = {
      id: `upgma_${nextId++}`,
      name: `Divergence: ${(minVal * 100).toFixed(3)}%`,
      isLeaf: false,
      height: newHeight,
      branchLength: 0,
      depthFromRoot: 0,
      pairwiseDist: minVal,
      indices: [...nodeA.indices, ...nodeB.indices],
      children: [nodeA, nodeB],
    };

    active.splice(v, 1);
    active.splice(u, 1);
    active.push(merged);
  }

  const root = active[0];
  const maxHeight = root.height || 0;
  const assignDepths = (node, parentDepth) => {
    if (node.id === root.id) {
      node.depthFromRoot = 0;
    } else {
      node.depthFromRoot = parentDepth + node.branchLength;
    }
    if (node.isLeaf) {
      node.height = maxHeight;
    }
    node.children.forEach((ch) => assignDepths(ch, node.depthFromRoot));
  };
  assignDepths(root, 0);
  return root;
}

/**
 * Saitou–Nei Neighbor-Joining followed by exact midpoint rooting on longest leaf-to-leaf path.
 * Negative branch lengths are clamped cleanly to 0.
 */
export function buildNeighborJoiningTree(sampleIds, matrix) {
  const nInit = Array.isArray(sampleIds) ? sampleIds.length : 0;
  if (nInit === 0) return null;
  if (nInit === 1) return buildUPGMATree(sampleIds, matrix);

  if (nInit === 2) {
    const d = Math.max(0, matrix[0]?.[1] ?? 0);
    const half = d / 2;
    return {
      id: 'nj_root',
      name: 'Midpoint Root',
      isLeaf: false,
      branchLength: 0,
      depthFromRoot: 0,
      height: half,
      pairwiseDist: d,
      children: [
        {
          id: 'leaf_0',
          name: sampleIds[0],
          isLeaf: true,
          branchLength: half,
          depthFromRoot: half,
          height: half,
          children: [],
        },
        {
          id: 'leaf_1',
          name: sampleIds[1],
          isLeaf: true,
          branchLength: half,
          depthFromRoot: half,
          height: half,
          children: [],
        },
      ],
    };
  }

  // Build undirected graph nodes and adjacency map
  let nextInternal = 1;
  const nodesMeta = new Map();
  const adj = new Map();

  const addNode = (id, name, isLeaf) => {
    nodesMeta.set(id, { id, name, isLeaf });
    adj.set(id, []);
  };
  const addEdge = (a, b, len) => {
    const w = Math.max(0, len);
    adj.get(a).push({ to: b, length: w });
    adj.get(b).push({ to: a, length: w });
  };
  const removeEdge = (a, b) => {
    adj.set(a, adj.get(a).filter((e) => e.to !== b));
    adj.set(b, adj.get(b).filter((e) => e.to !== a));
  };

  let active = [];
  for (let i = 0; i < nInit; i++) {
    const id = `leaf_${i}`;
    addNode(id, sampleIds[i], true);
    active.push(id);
  }

  // Mutable symmetric distance map keyed by active index
  let D = matrix.map((row) => row.map((v) => Math.max(0, v)));

  while (active.length > 2) {
    const r = active.length;
    const rowSums = Array(r).fill(0);
    for (let i = 0; i < r; i++) {
      for (let j = 0; j < r; j++) {
        if (i !== j) rowSums[i] += D[i][j];
      }
    }

    let minQ = Infinity;
    let bestI = 0;
    let bestJ = 1;
    for (let i = 0; i < r; i++) {
      for (let j = i + 1; j < r; j++) {
        const q = (r - 2) * D[i][j] - rowSums[i] - rowSums[j];
        if (q < minQ) {
          minQ = q;
          bestI = i;
          bestJ = j;
        }
      }
    }

    const dIJ = D[bestI][bestJ];
    const delta = (rowSums[bestI] - rowSums[bestJ]) / (2 * (r - 2));
    const lenI = Math.max(0, dIJ / 2 + delta);
    const lenJ = Math.max(0, dIJ - lenI);

    const newId = `nj_node_${nextInternal++}`;
    addNode(newId, newId, false);
    addEdge(newId, active[bestI], lenI);
    addEdge(newId, active[bestJ], lenJ);

    // Distance from new node to remaining nodes k
    const remainingIdx = [];
    for (let k = 0; k < r; k++) {
      if (k !== bestI && k !== bestJ) remainingIdx.push(k);
    }

    const m = remainingIdx.length;
    const newD = Array.from({ length: m + 1 }, () => Array(m + 1).fill(0));
    for (let a = 0; a < m; a++) {
      const ka = remainingIdx[a];
      for (let b = a + 1; b < m; b++) {
        const kb = remainingIdx[b];
        newD[a][b] = newD[b][a] = D[ka][kb];
      }
      const distNew = Math.max(0, (D[bestI][ka] + D[bestJ][ka] - dIJ) / 2);
      newD[a][m] = newD[m][a] = distNew;
    }

    active = remainingIdx.map((k) => active[k]);
    active.push(newId);
    D = newD;
  }

  // Join final two active nodes
  if (active.length === 2) {
    addEdge(active[0], active[1], D[0][1]);
  }

  // Midpoint rooting: find longest leaf-to-leaf path in undirected tree
  const leaves = [...nodesMeta.values()].filter((m) => m.isLeaf).map((m) => m.id);
  const findFarthest = (startId) => {
    const parent = new Map([[startId, null]]);
    const edgeToParent = new Map([[startId, 0]]);
    const dist = new Map([[startId, 0]]);
    const stack = [startId];
    while (stack.length > 0) {
      const cur = stack.pop();
      const curD = dist.get(cur);
      for (const edge of adj.get(cur) || []) {
        if (!dist.has(edge.to)) {
          dist.set(edge.to, curD + edge.length);
          parent.set(edge.to, cur);
          edgeToParent.set(edge.to, edge.length);
          stack.push(edge.to);
        }
      }
    }
    let farId = startId;
    let farD = -1;
    leaves.forEach((lid) => {
      const ld = dist.get(lid) ?? 0;
      if (ld > farD) {
        farD = ld;
        farId = lid;
      }
    });
    return { farId, farD, parent, edgeToParent };
  };

  const pass1 = findFarthest(leaves[0]);
  const pass2 = findFarthest(pass1.farId);

  // Reconstruct path from pass1.farId to pass2.farId
  const path = [];
  let cur = pass2.farId;
  while (cur !== null) {
    path.push({ id: cur, stepLen: pass2.edgeToParent.get(cur) || 0 });
    cur = pass2.parent.get(cur) ?? null;
  }
  path.reverse();

  const halfTotal = pass2.farD / 2;
  let accum = 0;
  let splitU = path[0].id;
  let splitV = path.length > 1 ? path[1].id : path[0].id;
  let lenFromU = 0;
  let lenFromV = 0;

  for (let i = 1; i < path.length; i++) {
    const seg = path[i].stepLen;
    if (accum + seg >= halfTotal || i === path.length - 1) {
      splitU = path[i - 1].id;
      splitV = path[i].id;
      lenFromU = Math.max(0, halfTotal - accum);
      lenFromV = Math.max(0, seg - lenFromU);
      break;
    }
    accum += seg;
  }

  addNode('nj_root', 'Midpoint Root', false);
  if (splitU !== splitV) {
    removeEdge(splitU, splitV);
    addEdge('nj_root', splitU, lenFromU);
    addEdge('nj_root', splitV, lenFromV);
  }

  // Orient tree from nj_root
  let internalSeq = 1;
  const orient = (nodeId, parentId, branchLen, depth) => {
    const meta = nodesMeta.get(nodeId);
    const childrenEdges = (adj.get(nodeId) || []).filter((e) => e.to !== parentId);
    const children = childrenEdges.map((e) => orient(e.to, nodeId, e.length, depth + e.length));
    const displayName =
      nodeId === 'nj_root'
        ? 'Midpoint Root'
        : meta.isLeaf
        ? meta.name
        : `Node ${internalSeq++}`;
    return {
      id: nodeId,
      name: displayName,
      isLeaf: meta.isLeaf,
      branchLength: branchLen,
      depthFromRoot: depth,
      children,
    };
  };

  const rooted = orient('nj_root', null, 0, 0);

  // Compute max leaf depth cleanly for consistent scale metadata
  const computeMaxDepth = (n) => {
    if (!n.children || n.children.length === 0) return n.depthFromRoot;
    return Math.max(...n.children.map(computeMaxDepth));
  };
  rooted.height = computeMaxDepth(rooted);
  rooted.pairwiseDist = pass2.farD;
  return rooted;
}
