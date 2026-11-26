import React, { useMemo, useRef, useState } from 'react';
import {
  AreaChart,
  Area,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Brush,
  ReferenceLine,
  ReferenceArea,
} from 'recharts';
import { AlertCircle, Download } from 'lucide-react';
import { getInternalDropoutIntervals } from '../utils/intervals';
import { serializeSvgElement, triggerBrowserDownload } from '../utils/exporters';
import { computeLogScaleConfig, formatDepthAxisTick } from '../utils/chartUtils';
import { formatDepth } from '../utils/reportSchema';

const CustomTooltip = ({ active, payload, label, minCoverageThreshold = 10 }) => {
  if (active && payload && payload.length) {
    const data = payload[0].payload;
    const depth = data.depth;

    let qualityLabel = 'FAIL';
    let badgeClass = 'bg-[#FCEEEB] text-[#A62C22] border-[#EBBDB8]';
    let depthClass = 'text-[#A62C22]';
    if (depth >= minCoverageThreshold) {
      qualityLabel = 'PASS';
      badgeClass = 'bg-[#ECF7F2] text-[#16654A] border-[#B5DFCE]';
      depthClass = 'text-[#16654A]';
    } else if (depth >= Math.max(1, Math.round(minCoverageThreshold * 0.25))) {
      qualityLabel = 'WARN';
      badgeClass = 'bg-[#FDF6E7] text-[#8C5308] border-[#EBD29B]';
      depthClass = 'text-[#8C5308]';
    }

    return (
      <div className="bg-white border border-[#CBD2D9] p-3 rounded-[6px] shadow-tooltip text-left select-text font-mono tabular-nums leading-relaxed max-w-xs">
        <div className="font-semibold text-[#111827] border-b border-[#E2E6EA] pb-1.5 mb-1.5 text-[11px] flex justify-between items-center gap-4">
          <span className="font-sans">Genomic Position</span>
          <span className={`text-[10px] px-1.5 py-0.5 rounded-[4px] border font-semibold ${badgeClass}`}>
            {qualityLabel}
          </span>
        </div>
        <div className="space-y-1 text-[11px] text-[#4B5563]">
          <p>
            Position: <span className="text-[#111827] font-semibold">{label ? label.toLocaleString() : 0} bp</span>
          </p>
          <p>
            Depth: <span className={`${depthClass} font-semibold`}>{formatDepth(depth, { decimals: 0 })}</span>
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
      <div className="flex flex-col items-center justify-center h-56 text-[#6B7280] glass-panel text-[13px]">
        <AlertCircle className="w-6 h-6 mb-2 text-[#A62C22]" />
        <p>No coverage profile available for this sample.</p>
      </div>
    );
  }

  const safeThreshold = Math.max(0, Number(minCoverageThreshold) || 0);
  const logConfig = computeLogScaleConfig(
    coverageData.map((d) => d.depth),
    safeThreshold
  );

  const chartData = coverageData.map((d) => ({
    ...d,
    displayDepth: selectedYScale === 'log' ? Math.log10(Math.max(1, d.depth)) : d.depth,
  }));

  const cutoffY = selectedYScale === 'log' ? Math.log10(Math.max(1, safeThreshold)) : safeThreshold;

  const handleExportSvg = () => {
    const svg = chartContainerRef.current?.querySelector('svg.recharts-surface');
    if (!svg) return;
    const svgText = serializeSvgElement(svg);
    triggerBrowserDownload(svgText, `viromics_${sampleId}_${selectedRefName}_coverage.svg`, 'image/svg+xml;charset=utf-8;');
  };

  return (
    <div className="p-5 glass-panel flex flex-col gap-4">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3">
        <div>
          <h3 className="text-[16px] font-semibold text-[#111827]">Genome Coverage Depth</h3>
          <p className="text-[12px] text-[#6B7280] mt-0.5">
            Drag the brush at the bottom to zoom into specific genes or amplicon zones
          </p>
        </div>

        <div className="flex items-center gap-2 shrink-0">
          <div className="flex bg-[#F1F3F5] border border-[#E2E6EA] rounded-[4px] p-0.5">
            <button
              onClick={() => setSelectedYScale('log')}
              className={`px-2.5 py-1 text-xs font-medium rounded-[3px] transition-colors ${
                selectedYScale === 'log'
                  ? 'bg-[#0F5C6E] text-white font-semibold'
                  : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              Log Scale
            </button>
            <button
              onClick={() => setSelectedYScale('linear')}
              className={`px-2.5 py-1 text-xs font-medium rounded-[3px] transition-colors ${
                selectedYScale === 'linear'
                  ? 'bg-[#0F5C6E] text-white font-semibold'
                  : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              Linear
            </button>
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

      {dropoutRegions.length > 0 && (
        <div className="bg-[#FCEEEB] border border-[#EBBDB8] rounded-[4px] px-3.5 py-2.5 flex items-start gap-2.5">
          <AlertCircle className="w-4 h-4 text-[#A62C22] shrink-0 mt-0.5" />
          <div className="flex-1 text-xs">
            <span className="font-semibold text-[#A62C22]">
              Internal Dropouts Detected ({dropoutRegions.length} {dropoutRegions.length === 1 ? 'region' : 'regions'})
            </span>
            <span className="text-[#4B5563] block mt-0.5 font-mono text-[11px] tabular-nums">
              Depth &lt; {minCoverageThreshold}x inside callable target intervals at:{' '}
              {dropoutRegions.slice(0, 4).map((r) => `${r.start.toLocaleString()}–${r.end.toLocaleString()}`).join(', ')}
              {dropoutRegions.length > 4 ? ` (+${dropoutRegions.length - 4} more)` : ''}. Masked as &apos;N&apos; in consensus FASTA.
            </span>
          </div>
        </div>
      )}

      <div ref={chartContainerRef} className="h-72 w-full mt-1 select-none">
        <ResponsiveContainer width="100%" height="100%">
          <AreaChart data={chartData} margin={{ top: 18, right: 16, left: 4, bottom: 0 }}>
            <defs>
              <linearGradient id="coverageGradient" x1="0" y1="0" x2="0" y2="1">
                <stop offset="5%" stopColor="#0F5C6E" stopOpacity={0.22} />
                <stop offset="95%" stopColor="#0F5C6E" stopOpacity={0.02} />
              </linearGradient>
            </defs>
            <CartesianGrid strokeDasharray="3 3" stroke="#E5E7EB" />
            {dropoutRegions.slice(0, 25).map((r, idx) => (
              <ReferenceArea
                key={`drop-${idx}`}
                x1={r.start}
                x2={r.end}
                fill="#C0392B"
                fillOpacity={0.1}
                ifOverflow="extendDomain"
              />
            ))}
            <XAxis
              dataKey="pos"
              tickFormatter={(tick) => (totalLength < 5000 ? `${tick}` : `${(tick / 1000).toFixed(1)}k`)}
              stroke="#CBD2D9"
              tick={{ fill: '#4B5563', fontSize: 11 }}
            />
            <YAxis
              width={48}
              domain={selectedYScale === 'log' ? logConfig.domain : [0, 'auto']}
              ticks={selectedYScale === 'log' ? logConfig.ticks : undefined}
              tickFormatter={(tick) => formatDepthAxisTick(tick, selectedYScale === 'log')}
              stroke="#CBD2D9"
              tick={{ fill: '#4B5563', fontSize: 11 }}
            />
            <Tooltip content={<CustomTooltip minCoverageThreshold={minCoverageThreshold} />} />
            <Area
              type="monotone"
              dataKey="displayDepth"
              stroke="#0F5C6E"
              strokeWidth={1.75}
              fillOpacity={1}
              fill="url(#coverageGradient)"
            />
            <ReferenceLine
              y={cutoffY}
              stroke="#C0392B"
              strokeDasharray="4 4"
              strokeWidth={1.25}
              label={{
                value: `Consensus Cutoff (${minCoverageThreshold}x)`,
                fill: '#A62C22',
                fontSize: 10,
                fontWeight: 600,
                position: 'insideTopRight',
              }}
            />
            <Brush
              dataKey="pos"
              height={28}
              stroke="#CBD2D9"
              fill="#F1F3F5"
              tickFormatter={(tick) => (totalLength < 5000 ? `${tick}` : `${(tick / 1000).toFixed(1)}k`)}
              travellerWidth={8}
            />
          </AreaChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
