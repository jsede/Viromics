import React, { useMemo, useRef, useState } from 'react';
import { AreaChart, Area, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Brush, ReferenceLine, ReferenceArea } from 'recharts';
import { AlertCircle, Download } from 'lucide-react';
import { getInternalDropoutIntervals } from '../utils/intervals';
import { serializeSvgElement, triggerBrowserDownload } from '../utils/exporters';

const CustomTooltip = ({ active, payload, label, minCoverageThreshold = 10 }) => {
  if (active && payload && payload.length) {
    const data = payload[0].payload;
    const depth = data.depth;

    let qualityLabel = 'FAIL';
    let qualityColor = 'text-rose-400';
    if (depth >= minCoverageThreshold) {
      qualityLabel = 'PASS';
      qualityColor = 'text-emerald-400';
    } else if (depth >= Math.max(1, Math.round(minCoverageThreshold * 0.25))) {
      qualityLabel = 'WARN';
      qualityColor = 'text-amber-400';
    }

    return (
      <div className="bg-slate-900 border border-slate-700/80 p-3 rounded-xl shadow-xl text-left select-text font-mono leading-relaxed max-w-xs">
        <div className="font-bold text-white border-b border-slate-800 pb-1.5 mb-1.5 text-[10px] tracking-wider flex justify-between items-center gap-4">
          <span>Genomic Position</span>
          <span className={`text-[9px] px-1.5 py-0.5 rounded bg-slate-800 font-bold ${qualityColor}`}>
            {qualityLabel}
          </span>
        </div>
        <div className="space-y-1 text-[11px]">
          <p className="text-slate-300">
            Position: <span className="text-white font-bold">{label ? label.toLocaleString() : 0} bp</span>
          </p>
          <p className="text-slate-300">
            Depth: <span className={`${qualityColor} font-bold`}>{depth.toLocaleString()}x</span>
          </p>
        </div>
      </div>
    );
  }
  return null;
};

export default function CoverageChart({
  coverageData,
  minCoverageThreshold = 10,
  lowCoverageIntervals = null,
  callableTargetIntervals = null,
  sampleId = 'sample',
  selectedRefName = 'contig',
}) {
  const [selectedYScale, setSelectedYScale] = useState('log');
  const chartContainerRef = useRef(null);

  const totalLength = coverageData && coverageData.length > 0 ? coverageData[coverageData.length - 1].pos : 0;

  const dropoutRegions = useMemo(
    () =>
      getInternalDropoutIntervals(
        lowCoverageIntervals,
        callableTargetIntervals,
        coverageData,
        minCoverageThreshold
      ),
    [lowCoverageIntervals, callableTargetIntervals, coverageData, minCoverageThreshold]
  );

  if (!coverageData || coverageData.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center h-64 text-slate-500 glass-panel rounded-xl">
        <AlertCircle className="w-8 h-8 mb-2 text-red-600 animate-pulse" />
        <p>No coverage profile available for this sample.</p>
      </div>
    );
  }

  const chartData = coverageData.map((d) => ({
    ...d,
    displayDepth: selectedYScale === 'log' ? Math.log10(Math.max(1, d.depth)) : d.depth,
  }));

  const safeThreshold = Math.max(0, Number(minCoverageThreshold) || 0);
  const cutoffY = selectedYScale === 'log' ? Math.log10(Math.max(1, safeThreshold)) : safeThreshold;

  const formatYAxis = (tick) => {
    if (selectedYScale === 'log') {
      const val = Math.round(Math.pow(10, tick));
      return val >= 1000 ? `${val / 1000}k` : val;
    }
    return tick >= 1000 ? `${tick / 1000}k` : tick;
  };

  const handleExportSvg = () => {
    const svg = chartContainerRef.current?.querySelector('svg.recharts-surface');
    if (!svg) return;
    const svgText = serializeSvgElement(svg);
    triggerBrowserDownload(svgText, `viromics_${sampleId}_${selectedRefName}_coverage.svg`, 'image/svg+xml;charset=utf-8;');
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-4">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold text-slate-900">Genome Coverage Depth</h3>
          <p className="text-xs text-slate-500">
            Drag the brush at the bottom to zoom into specific genes or amplicon zones
          </p>
        </div>

        <div className="flex items-center gap-2">
          <div className="flex bg-slate-100 border border-slate-200 rounded-lg p-0.5">
            <button
              onClick={() => setSelectedYScale('log')}
              className={`px-3 py-1 text-xs font-medium rounded-md transition-all ${
                selectedYScale === 'log'
                  ? 'bg-red-600 text-white font-semibold'
                  : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              Log Scale
            </button>
            <button
              onClick={() => setSelectedYScale('linear')}
              className={`px-3 py-1 text-xs font-medium rounded-md transition-all ${
                selectedYScale === 'linear'
                  ? 'bg-red-600 text-white font-semibold'
                  : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              Linear
            </button>
          </div>
          <button
            onClick={handleExportSvg}
            className="flex items-center gap-1.5 px-3 py-1.5 bg-white border border-slate-300 text-slate-700 text-xs font-semibold rounded-lg hover:bg-slate-50 transition-colors"
          >
            <Download className="w-3.5 h-3.5 text-red-600" />
            Export SVG
          </button>
        </div>
      </div>

      {dropoutRegions.length > 0 && (
        <div className="bg-red-50 border border-red-200 rounded-xl p-3 flex items-start gap-3">
          <AlertCircle className="w-5 h-5 text-red-600 shrink-0 mt-0.5" />
          <div className="flex-1">
            <span className="text-xs font-semibold text-red-700 block">
              Internal Dropouts Detected ({dropoutRegions.length} {dropoutRegions.length === 1 ? 'region' : 'regions'})
            </span>
            <span className="text-[11px] text-slate-700 block mt-1">
              Depth fell below the {minCoverageThreshold}x consensus threshold inside callable target intervals at:{' '}
              {dropoutRegions.slice(0, 4).map((r) => `${r.start.toLocaleString()}–${r.end.toLocaleString()}`).join(', ')}
              {dropoutRegions.length > 4 ? ` (+${dropoutRegions.length - 4} more)` : ''}.
              These positions are masked as &apos;N&apos; in the consensus FASTA.
            </span>
          </div>
        </div>
      )}

      <div ref={chartContainerRef} className="h-72 w-full mt-2 select-none">
        <ResponsiveContainer width="100%" height="100%">
          <AreaChart data={chartData} margin={{ top: 10, right: 10, left: -20, bottom: 0 }}>
            <defs>
              <linearGradient id="coverageGradient" x1="0" y1="0" x2="0" y2="1">
                <stop offset="5%" stopColor="#DC2626" stopOpacity={0.35} />
                <stop offset="95%" stopColor="#DC2626" stopOpacity={0.02} />
              </linearGradient>
            </defs>
            <CartesianGrid strokeDasharray="3 3" stroke="#E2E8F0" />
            {dropoutRegions.slice(0, 25).map((r, idx) => (
              <ReferenceArea
                key={`drop-${idx}`}
                x1={r.start}
                x2={r.end}
                fill="#EF4444"
                fillOpacity={0.12}
                ifOverflow="extendDomain"
              />
            ))}
            <XAxis
              dataKey="pos"
              tickFormatter={(tick) => (totalLength < 5000 ? `${tick}` : `${(tick / 1000).toFixed(1)}k`)}
              stroke="#CBD5E1"
              tick={{ fill: '#475569', fontSize: 10 }}
            />
            <YAxis
              tickFormatter={formatYAxis}
              stroke="#CBD5E1"
              tick={{ fill: '#475569', fontSize: 10 }}
            />
            <Tooltip content={<CustomTooltip minCoverageThreshold={minCoverageThreshold} />} />
            <ReferenceLine
              y={cutoffY}
              stroke="#DC2626"
              strokeDasharray="4 4"
              label={{
                value: `Consensus Cutoff (${minCoverageThreshold}x)`,
                fill: '#DC2626',
                fontSize: 10,
                position: 'top',
              }}
            />
            <Area
              type="monotone"
              dataKey="displayDepth"
              stroke="#DC2626"
              strokeWidth={1.75}
              fillOpacity={1}
              fill="url(#coverageGradient)"
            />
            <Brush
              dataKey="pos"
              height={30}
              stroke="#94A3B8"
              fill="#F1F5F9"
              tickFormatter={(tick) => (totalLength < 5000 ? `${tick}` : `${(tick / 1000).toFixed(1)}k`)}
              travellerWidth={8}
            />
          </AreaChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
