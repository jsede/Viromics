import React, { useMemo, useRef, useState } from 'react';
import {
  BarChart,
  Bar,
  Cell,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  ReferenceLine,
} from 'recharts';
import { Download } from 'lucide-react';
import { serializeSvgElement, triggerBrowserDownload } from '../utils/exporters';
import { computeLogScaleConfig, formatDepthAxisTick } from '../utils/chartUtils';
import { formatDepth } from '../utils/reportSchema';

const CustomAmpliconTooltip = ({ active, payload, minCoverageThreshold }) => {
  if (!active || !payload || payload.length === 0) return null;
  const amp = payload[0].payload;
  const isPass = amp.status === 'pass';

  return (
    <div className="bg-white border border-[#CBD2D9] p-3 rounded-[6px] shadow-tooltip text-left font-mono tabular-nums text-[11px] text-[#4B5563] max-w-xs">
      <div className="font-semibold text-[#111827] border-b border-[#E2E6EA] pb-1.5 mb-1.5 flex justify-between items-center gap-3">
        <span>{amp.amplicon_id}</span>
        <span
          className={`px-1.5 py-0.5 rounded-[4px] border text-[10px] font-semibold uppercase ${
            isPass
              ? 'bg-[#ECF7F2] text-[#16654A] border-[#B5DFCE]'
              : 'bg-[#FCEEEB] text-[#A62C22] border-[#EBBDB8]'
          }`}
        >
          {amp.status}
        </span>
      </div>
      <div className="space-y-0.5">
        <p>
          Insert Span: <strong className="text-[#111827]">{amp.start.toLocaleString()}–{amp.end.toLocaleString()} bp</strong>
        </p>
        <p>
          Mean Depth:{' '}
          <strong className={isPass ? 'text-[#16654A]' : 'text-[#A62C22]'}>
            {formatDepth(amp.mean_depth)}
          </strong>
        </p>
        <p>
          Median Depth: <strong className="text-[#111827]">{formatDepth(amp.median_depth)}</strong>
        </p>
        <p>
          Min Depth: <strong className="text-[#111827]">{formatDepth(amp.min_depth, { decimals: 0 })}</strong>
        </p>
        <p>
          Bases &ge; {minCoverageThreshold}x:{' '}
          <strong className={amp.percent_bases_above_threshold >= 90 ? 'text-[#16654A]' : 'text-[#8C5308]'}>
            {amp.percent_bases_above_threshold.toFixed(1)}%
          </strong>
        </p>
      </div>
    </div>
  );
};

function ThresholdBadgeLabel({ viewBox, value }) {
  if (!viewBox) return null;
  const textWidth = Math.max(46, String(value).length * 6.5 + 10);
  const pillX = viewBox.x + viewBox.width - textWidth - 4;
  const pillY = Math.max(4, viewBox.y - 19);
  return (
    <g>
      <rect
        x={pillX}
        y={pillY}
        width={textWidth}
        height={16}
        rx="3"
        fill="#FFFFFF"
        fillOpacity="0.95"
        stroke="#EBBDB8"
        strokeWidth="1"
      />
      <text
        x={pillX + textWidth / 2}
        y={pillY + 11.5}
        fill="#A62C22"
        fontSize="10"
        fontWeight="600"
        textAnchor="middle"
      >
        {value}
      </text>
    </g>
  );
}

export default function AmpliconChart({
  amplicons = [],
  selectedRefName = '',
  sampleId = '',
  minCoverageThreshold = 10,
}) {
  const [yScale, setYScale] = useState('log');
  const [statusFilter, setStatusFilter] = useState('ALL');
  const containerRef = useRef(null);

  const refAmplicons = useMemo(() => {
    if (!Array.isArray(amplicons)) return [];
    const matching = selectedRefName
      ? amplicons.filter((a) => a.reference === selectedRefName)
      : amplicons;
    return [...matching].sort((a, b) => a.start - b.start);
  }, [amplicons, selectedRefName]);

  const filteredAmplicons = useMemo(() => {
    if (statusFilter === 'ALL') return refAmplicons;
    return refAmplicons.filter((a) => a.status === statusFilter);
  }, [refAmplicons, statusFilter]);

  if (refAmplicons.length === 0) {
    return null;
  }

  const passingCount = refAmplicons.filter((a) => a.status === 'pass').length;
  const failingCount = refAmplicons.length - passingCount;
  const chartData = filteredAmplicons.map((a) => ({
    ...a,
    displayDepth: yScale === 'log' ? Math.log10(Math.max(1, a.mean_depth)) : a.mean_depth,
  }));

  const logConfig = computeLogScaleConfig(
    refAmplicons.map((a) => a.mean_depth),
    minCoverageThreshold
  );
  const safeRefCutoff =
    yScale === 'log' ? Math.log10(Math.max(1, minCoverageThreshold)) : Math.max(0, minCoverageThreshold);

  const handleExportSvg = () => {
    const svg = containerRef.current?.querySelector('svg.recharts-surface');
    if (!svg) return;
    const svgText = serializeSvgElement(svg);
    triggerBrowserDownload(svgText, `viromics_${sampleId}_${selectedRefName || 'contig'}_amplicons.svg`, 'image/svg+xml;charset=utf-8;');
  };

  return (
    <div className="p-5 glass-panel flex flex-col gap-4">
      <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-3">
        <div className="min-w-0">
          <h3 className="text-[16px] font-semibold text-[#111827]">Amplicon Panel QC</h3>
          <p className="text-[12px] text-[#6B7280] mt-0.5 tabular-nums">
            Insert coverage ({refAmplicons.length} amplicons:{' '}
            <strong className="text-[#16654A] font-semibold">{passingCount} pass</strong>,{' '}
            <strong className={failingCount > 0 ? 'text-[#A62C22] font-semibold' : 'text-[#4B5563] font-semibold'}>
              {failingCount} fail
            </strong>
            ; pass &ge;90% bases &ge;{minCoverageThreshold}x &amp; mean &ge;{minCoverageThreshold}x)
          </p>
        </div>

        <div className="flex items-center gap-2 shrink-0 flex-nowrap overflow-x-auto">
          <div className="flex bg-[#F1F3F5] border border-[#E2E6EA] rounded-[4px] p-0.5 text-xs">
            {[
              { id: 'ALL', label: `All (${refAmplicons.length})` },
              { id: 'pass', label: `Pass (${passingCount})` },
              { id: 'fail', label: `Fail (${failingCount})` },
            ].map((opt) => (
              <button
                key={opt.id}
                onClick={() => setStatusFilter(opt.id)}
                className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                  statusFilter === opt.id
                    ? 'bg-[#0F5C6E] text-white font-semibold'
                    : 'text-[#4B5563] hover:text-[#111827]'
                }`}
              >
                {opt.label}
              </button>
            ))}
          </div>

          <div className="flex bg-[#F1F3F5] border border-[#E2E6EA] rounded-[4px] p-0.5 text-xs">
            <button
              onClick={() => setYScale('log')}
              className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                yScale === 'log' ? 'bg-[#0F5C6E] text-white font-semibold' : 'text-[#4B5563] hover:text-[#111827]'
              }`}
            >
              Log Scale
            </button>
            <button
              onClick={() => setYScale('linear')}
              className={`px-2.5 py-1 rounded-[3px] font-medium transition-colors whitespace-nowrap ${
                yScale === 'linear' ? 'bg-[#0F5C6E] text-white font-semibold' : 'text-[#4B5563] hover:text-[#111827]'
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

      <div ref={containerRef} className="h-64 w-full select-none">
        <ResponsiveContainer width="100%" height="100%">
          <BarChart data={chartData} margin={{ top: 20, right: 16, left: 4, bottom: 20 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#E5E7EB" />
            <XAxis
              dataKey="amplicon_id"
              interval={0}
              angle={chartData.length > 10 ? -35 : 0}
              textAnchor={chartData.length > 10 ? 'end' : 'middle'}
              height={chartData.length > 10 ? 55 : 28}
              stroke="#CBD2D9"
              tick={{ fill: '#4B5563', fontSize: 10 }}
            />
            <YAxis
              width={48}
              domain={yScale === 'log' ? logConfig.domain : [0, 'auto']}
              ticks={yScale === 'log' ? logConfig.ticks : undefined}
              tickFormatter={(tick) => formatDepthAxisTick(tick, yScale === 'log')}
              stroke="#CBD2D9"
              tick={{ fill: '#4B5563', fontSize: 11 }}
            />
            <Tooltip content={<CustomAmpliconTooltip minCoverageThreshold={minCoverageThreshold} />} />
            <Bar dataKey="displayDepth" radius={[3, 3, 0, 0]}>
              {chartData.map((entry, idx) => (
                <Cell key={`amp-${idx}`} fill={entry.status === 'pass' ? '#2A8C6A' : '#C0392B'} />
              ))}
            </Bar>
            <ReferenceLine
              y={safeRefCutoff}
              stroke="#C0392B"
              strokeDasharray="4 4"
              strokeWidth={1.25}
              ifOverflow="extendDomain"
              label={<ThresholdBadgeLabel value={`${minCoverageThreshold}x Min`} />}
            />
          </BarChart>
        </ResponsiveContainer>
      </div>

      <div className="overflow-x-auto max-h-56 overflow-y-auto border border-[#E2E6EA] rounded-[6px]">
        <table className="w-full text-left text-xs border-collapse">
          <thead className="sticky top-0 bg-[#F1F3F5] border-b border-[#E2E6EA] text-[#4B5563] text-[11px] uppercase tracking-[0.04em]">
            <tr>
              <th className="px-3 py-2 font-semibold">Amplicon</th>
              <th className="px-3 py-2 font-semibold">Insert Coordinates</th>
              <th className="px-3 py-2 font-semibold text-right">Mean Depth</th>
              <th className="px-3 py-2 font-semibold text-right">Median</th>
              <th className="px-3 py-2 font-semibold text-right">Min Depth</th>
              <th className="px-3 py-2 font-semibold text-right">% &ge; {minCoverageThreshold}x</th>
              <th className="px-3 py-2 font-semibold text-center">Status</th>
            </tr>
          </thead>
          <tbody>
            {filteredAmplicons.map((a) => (
              <tr key={a.amplicon_id} className="border-b border-[#E2E6EA] last:border-none bg-white hover:bg-[#F8F9FA]">
                <td className="px-3 py-2 font-mono font-semibold text-[#111827]">{a.amplicon_id}</td>
                <td className="px-3 py-2 font-mono text-[#4B5563] whitespace-nowrap">
                  {a.start.toLocaleString()}–{a.end.toLocaleString()} bp
                </td>
                <td className="px-3 py-2 font-mono text-right font-semibold text-[#111827] whitespace-nowrap">
                  {formatDepth(a.mean_depth)}
                </td>
                <td className="px-3 py-2 font-mono text-right text-[#4B5563] whitespace-nowrap">
                  {formatDepth(a.median_depth)}
                </td>
                <td className="px-3 py-2 font-mono text-right text-[#4B5563] whitespace-nowrap">
                  {formatDepth(a.min_depth, { decimals: 0 })}
                </td>
                <td className="px-3 py-2 font-mono text-right font-semibold text-[#111827] whitespace-nowrap">
                  {a.percent_bases_above_threshold.toFixed(1)}%
                </td>
                <td className="px-3 py-2 text-center whitespace-nowrap">
                  <span
                    className={`inline-block px-2 py-0.5 rounded-[4px] border text-[10px] font-semibold uppercase ${
                      a.status === 'pass'
                        ? 'bg-[#ECF7F2] text-[#16654A] border-[#B5DFCE]'
                        : 'bg-[#FCEEEB] text-[#A62C22] border-[#EBBDB8]'
                    }`}
                  >
                    {a.status}
                  </span>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}
