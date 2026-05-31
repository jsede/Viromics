import React, { useState } from 'react';
import { AreaChart, Area, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Brush, ReferenceLine } from 'recharts';
import { ZoomIn, ZoomOut, AlertCircle } from 'lucide-react';

export default function CoverageChart({ coverageData, minCoverageThreshold = 10 }) {
  const [zoomRange, setZoomRange] = useState({ start: 0, end: coverageData.length - 1 });
  const [selectedYScale, setSelectedYScale] = useState('log'); // 'linear' or 'log'

  if (!coverageData || coverageData.length === 0) {
    return (
      <div className="flex flex-col items-center justify-center h-64 text-gray-400 glass-panel rounded-xl">
        <AlertCircle className="w-8 h-8 mb-2 text-brand-teal animate-pulse" />
        <p>No coverage profile available for this sample.</p>
      </div>
    );
  }

  // Find dropout regions (consecutive positions where depth < threshold)
  // We ignore dropouts that occur entirely at the extreme terminal ends (first/last 50 bp)
  // as these are expected primer alignment boundaries rather than true amplicon failed reactions.
  const dropoutRegions = [];
  let currentDropout = null;
  const terminalBuffer = 50;
  const totalLength = coverageData.length > 0 ? coverageData[coverageData.length - 1].pos : 0;

  coverageData.forEach((d, idx) => {
    if (d.depth < minCoverageThreshold) {
      if (currentDropout === null) {
        currentDropout = { start: d.pos, end: d.pos };
      } else {
        currentDropout.end = d.pos;
      }
    } else {
      if (currentDropout !== null) {
        const isTerminal = currentDropout.end <= terminalBuffer || currentDropout.start >= (totalLength - terminalBuffer);
        if (!isTerminal) {
          dropoutRegions.push(currentDropout);
        }
        currentDropout = null;
      }
    }
  });
  if (currentDropout !== null) {
    const isTerminal = currentDropout.end <= terminalBuffer || currentDropout.start >= (totalLength - terminalBuffer);
    if (!isTerminal) {
      dropoutRegions.push(currentDropout);
    }
  }

  // Pre-process scale (avoid log(0))
  const chartData = coverageData.map(d => ({
    ...d,
    displayDepth: selectedYScale === 'log' ? Math.log10(Math.max(1, d.depth)) : d.depth,
  }));

  // Tooltip formatter
  const formatTooltip = (value, name, props) => {
    const rawDepth = props.payload.depth;
    return [`${rawDepth.toLocaleString()}x`, 'Depth'];
  };

  const formatYAxis = (tick) => {
    if (selectedYScale === 'log') {
      const val = Math.round(Math.pow(10, tick));
      return val >= 1000 ? `${val/1000}k` : val;
    }
    return tick >= 1000 ? `${tick/1000}k` : tick;
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-4">
      {/* Chart Header */}
      <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold font-outfit text-white">Genome Coverage Depth</h3>
          <p className="text-xs text-gray-400">
            Drag the brush at the bottom to zoom into specific genes or amplicon zones
          </p>
        </div>
        
        <div className="flex items-center gap-3">
          <div className="flex bg-dark-bg border border-dark-border rounded-lg p-0.5">
            <button
              onClick={() => setSelectedYScale('log')}
              className={`px-3 py-1 text-xs font-medium rounded-md transition-all ${
                selectedYScale === 'log'
                  ? 'bg-brand-teal text-dark-bg font-semibold'
                  : 'text-gray-400 hover:text-white'
              }`}
            >
              Log Scale
            </button>
            <button
              onClick={() => setSelectedYScale('linear')}
              className={`px-3 py-1 text-xs font-medium rounded-md transition-all ${
                selectedYScale === 'linear'
                  ? 'bg-brand-teal text-dark-bg font-semibold'
                  : 'text-gray-400 hover:text-white'
              }`}
            >
              Linear
            </button>
          </div>
        </div>
      </div>

      {/* Coverage Warnings / Dropout Alerts */}
      {dropoutRegions.length > 0 && (
        <div className="bg-red-50 border border-red-200 rounded-xl p-3 flex items-start gap-3">
          <AlertCircle className="w-5 h-5 text-brand-rose shrink-0 mt-0.5" />
          <div className="flex-1">
            <span className="text-xs font-semibold text-brand-rose block font-outfit">
              Amplicon Dropouts Detected! ({dropoutRegions.length} regions)
            </span>
            <span className="text-[11px] text-gray-600 block mt-1">
              Coverage fell below the {minCoverageThreshold}x consensus threshold in regions: {' '}
              {dropoutRegions.slice(0, 3).map(r => `${r.start.toLocaleString()}-${r.end.toLocaleString()}`).join(', ')}
              {dropoutRegions.length > 3 ? ` (+${dropoutRegions.length - 3} more)` : ''}.
              These bases will be masked as 'N' in the consensus FASTA.
            </span>
          </div>
        </div>
      )}

      {/* Recharts Coverage Area Chart */}
      <div className="h-72 w-full mt-2 select-none">
        <ResponsiveContainer width="100%" height="100%">
          <AreaChart
            data={chartData}
            margin={{ top: 10, right: 10, left: -20, bottom: 0 }}
          >
            <defs>
              <linearGradient id="coverageGradient" x1="0" y1="0" x2="0" y2="1">
                <stop offset="5%" stopColor="#14b8a6" stopOpacity={0.4}/>
                <stop offset="95%" stopColor="#14b8a6" stopOpacity={0.0}/>
              </linearGradient>
            </defs>
            <CartesianGrid strokeDasharray="3 3" stroke="rgba(36, 47, 77, 0.3)" />
            <XAxis 
              dataKey="pos" 
              tickFormatter={(tick) => totalLength < 5000 ? `${tick}` : `${(tick/1000).toFixed(1)}k`}
              stroke="rgba(156, 163, 175, 0.5)"
              tick={{ fill: '#9ca3af', fontSize: 10 }}
            />
            <YAxis 
              tickFormatter={formatYAxis}
              stroke="rgba(156, 163, 175, 0.5)"
              tick={{ fill: '#9ca3af', fontSize: 10 }}
            />
            <Tooltip
              contentStyle={{
                backgroundColor: 'rgba(22, 29, 48, 0.95)',
                borderColor: 'rgba(36, 47, 77, 0.8)',
                borderRadius: '12px',
                color: '#fff',
                fontFamily: 'Fira Code',
                fontSize: '11px'
              }}
              formatter={formatTooltip}
              labelFormatter={(label) => `Genomic Position: ${label.toLocaleString()} bp`}
            />
            <ReferenceLine 
              y={selectedYScale === 'log' ? Math.log10(minCoverageThreshold) : minCoverageThreshold} 
              stroke="#f43f5e" 
              strokeDasharray="4 4"
              label={{
                value: `Consensus Cutoff (${minCoverageThreshold}x)`,
                fill: '#f43f5e',
                fontSize: 10,
                position: 'top'
              }}
            />
            <Area
              type="monotone"
              dataKey="displayDepth"
              stroke="#14b8a6"
              strokeWidth={1.5}
              fillOpacity={1}
              fill="url(#coverageGradient)"
            />
            <Brush 
              dataKey="pos" 
              height={30} 
              stroke="#242f4d"
              fill="#161d30"
              tickFormatter={(tick) => totalLength < 5000 ? `${tick}` : `${(tick/1000).toFixed(1)}k`}
              travellerWidth={8}
            />
          </AreaChart>
        </ResponsiveContainer>
      </div>
    </div>
  );
}
