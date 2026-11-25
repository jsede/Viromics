import React, { useMemo, useRef, useState } from 'react';
import { BarChart, Bar, Cell, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ReferenceLine } from 'recharts';
import { Activity, Download } from 'lucide-react';
import { serializeSvgElement, triggerBrowserDownload } from '../utils/exporters';

const CustomAmpliconTooltip = ({ active, payload, minCoverageThreshold }) => {
  if (!active || !payload || payload.length === 0) return null;
  const amp = payload[0].payload;
  const isPass = amp.status === 'pass';

  return (
    <div className="bg-slate-900 border border-slate-700 p-3 rounded-xl shadow-xl text-left font-mono text-[11px] text-slate-200 max-w-xs">
      <div className="font-bold text-white border-b border-slate-800 pb-1.5 mb-1.5 flex justify-between items-center gap-3">
        <span>{amp.amplicon_id}</span>
        <span className={`px-1.5 py-0.5 rounded text-[9px] font-bold uppercase ${isPass ? 'bg-emerald-950 text-emerald-400' : 'bg-red-950 text-red-400'}`}>
          {amp.status}
        </span>
      </div>
      <div className="space-y-0.5">
        <p>Insert Span: <strong className="text-white">{amp.start.toLocaleString()}–{amp.end.toLocaleString()} bp</strong></p>
        <p>Mean Depth: <strong className={isPass ? 'text-emerald-400' : 'text-red-400'}>{amp.mean_depth.toLocaleString()}x</strong></p>
        <p>Median Depth: <strong className="text-white">{amp.median_depth.toLocaleString()}x</strong></p>
        <p>Min Depth: <strong className="text-white">{amp.min_depth.toLocaleString()}x</strong></p>
        <p>Bases &ge; {minCoverageThreshold}x: <strong className={amp.percent_bases_above_threshold >= 90 ? 'text-emerald-400' : 'text-amber-400'}>{amp.percent_bases_above_threshold.toFixed(1)}%</strong></p>
      </div>
    </div>
  );
};

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

  const safeRefCutoff = yScale === 'log' ? Math.log10(Math.max(1, minCoverageThreshold)) : Math.max(0, minCoverageThreshold);

  const formatYAxis = (tick) => {
    const val = yScale === 'log' ? Math.round(Math.pow(10, tick)) : tick;
    return val >= 1000 ? `${(val / 1000).toFixed(val % 1000 === 0 ? 0 : 1)}k` : String(val);
  };

  const handleExportSvg = () => {
    const svg = containerRef.current?.querySelector('svg.recharts-surface');
    if (!svg) return;
    const svgText = serializeSvgElement(svg);
    triggerBrowserDownload(svgText, `viromics_${sampleId}_${selectedRefName || 'contig'}_amplicons.svg`, 'image/svg+xml;charset=utf-8;');
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-5">
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <div className="flex items-center gap-2">
            <Activity className="w-5 h-5 text-red-600" />
            <h3 className="text-lg font-semibold text-slate-900">Amplicon Panel QC</h3>
          </div>
          <p className="text-xs text-slate-500 mt-0.5">
            Primer-trimmed insert coverage ({refAmplicons.length} amplicons: {' '}
            <strong className="text-emerald-700">{passingCount} pass</strong>, {' '}
            <strong className={failingCount > 0 ? 'text-red-600' : 'text-slate-600'}>{failingCount} fail</strong>; pass requires &ge;90% bases &ge;{minCoverageThreshold}x &amp; mean depth &ge;{minCoverageThreshold}x)
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-2">
          <div className="flex bg-slate-100 border border-slate-200 rounded-lg p-0.5 text-xs">
            {[
              { id: 'ALL', label: `All (${refAmplicons.length})` },
              { id: 'pass', label: `Pass (${passingCount})` },
              { id: 'fail', label: `Fail (${failingCount})` },
            ].map((opt) => (
              <button
                key={opt.id}
                onClick={() => setStatusFilter(opt.id)}
                className={`px-2.5 py-1 rounded-md font-medium transition-colors ${
                  statusFilter === opt.id
                    ? 'bg-red-600 text-white font-semibold'
                    : 'text-slate-600 hover:text-slate-900'
                }`}
              >
                {opt.label}
              </button>
            ))}
          </div>

          <div className="flex bg-slate-100 border border-slate-200 rounded-lg p-0.5 text-xs">
            <button
              onClick={() => setYScale('log')}
              className={`px-3 py-1 rounded-md font-medium transition-colors ${
                yScale === 'log' ? 'bg-red-600 text-white font-semibold' : 'text-slate-600 hover:text-slate-900'
              }`}
            >
              Log Scale
            </button>
            <button
              onClick={() => setYScale('linear')}
              className={`px-3 py-1 rounded-md font-medium transition-colors ${
                yScale === 'linear' ? 'bg-red-600 text-white font-semibold' : 'text-slate-600 hover:text-slate-900'
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

      <div ref={containerRef} className="h-64 w-full select-none">
        <ResponsiveContainer width="100%" height="100%">
          <BarChart data={chartData} margin={{ top: 12, right: 12, left: -12, bottom: 20 }}>
            <CartesianGrid strokeDasharray="3 3" stroke="#E2E8F0" />
            <XAxis
              dataKey="amplicon_id"
              interval={0}
              angle={chartData.length > 10 ? -35 : 0}
              textAnchor={chartData.length > 10 ? 'end' : 'middle'}
              height={chartData.length > 10 ? 55 : 28}
              tick={{ fill: '#475569', fontSize: 10 }}
            />
            <YAxis tickFormatter={formatYAxis} tick={{ fill: '#475569', fontSize: 10 }} />
            <Tooltip content={<CustomAmpliconTooltip minCoverageThreshold={minCoverageThreshold} />} />
            <ReferenceLine
              y={safeRefCutoff}
              stroke="#DC2626"
              strokeDasharray="4 4"
              label={{
                value: `${minCoverageThreshold}x Min`,
                fill: '#DC2626',
                fontSize: 10,
                position: 'top',
              }}
            />
            <Bar dataKey="displayDepth" radius={[4, 4, 0, 0]}>
              {chartData.map((entry, idx) => (
                <Cell key={`amp-${idx}`} fill={entry.status === 'pass' ? '#10B981' : '#EF4444'} />
              ))}
            </Bar>
          </BarChart>
        </ResponsiveContainer>
      </div>

      <div className="overflow-x-auto max-h-56 overflow-y-auto border border-slate-200 rounded-xl">
        <table className="w-full text-left text-xs border-collapse">
          <thead className="sticky top-0 bg-slate-100 border-b border-slate-200 text-slate-700">
            <tr>
              <th className="p-2.5 font-semibold">Amplicon</th>
              <th className="p-2.5 font-semibold">Insert Coordinates</th>
              <th className="p-2.5 font-semibold text-right">Mean Depth</th>
              <th className="p-2.5 font-semibold text-right">Median</th>
              <th className="p-2.5 font-semibold text-right">Min Depth</th>
              <th className="p-2.5 font-semibold text-right">% &ge; {minCoverageThreshold}x</th>
              <th className="p-2.5 font-semibold text-center">Status</th>
            </tr>
          </thead>
          <tbody>
            {filteredAmplicons.map((a) => (
              <tr key={a.amplicon_id} className="border-b border-slate-100 bg-white hover:bg-slate-50">
                <td className="p-2.5 font-mono font-semibold text-slate-900">{a.amplicon_id}</td>
                <td className="p-2.5 font-mono text-slate-600">{a.start.toLocaleString()}–{a.end.toLocaleString()} bp</td>
                <td className="p-2.5 font-mono text-right font-semibold text-slate-800">{a.mean_depth.toLocaleString()}x</td>
                <td className="p-2.5 font-mono text-right text-slate-600">{a.median_depth.toLocaleString()}x</td>
                <td className="p-2.5 font-mono text-right text-slate-600">{a.min_depth.toLocaleString()}x</td>
                <td className="p-2.5 font-mono text-right font-semibold text-slate-800">{a.percent_bases_above_threshold.toFixed(1)}%</td>
                <td className="p-2.5 text-center">
                  <span className={`inline-block px-2 py-0.5 rounded-full text-[10px] font-bold uppercase ${
                    a.status === 'pass'
                      ? 'bg-emerald-100 text-emerald-700'
                      : 'bg-red-100 text-red-700'
                  }`}>
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
