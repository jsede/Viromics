import React, { useState, useEffect } from 'react';
import {
  Dna,
  UploadCloud,
  Grid,
  Sliders,
  Activity,
  BarChart3,
  GitBranch,
  AlertTriangle,
  CheckCircle2,
  XCircle,
  FileText,
  Layers,
  ShieldCheck,
  Zap,
  Printer,
  ArrowUpDown,
  Download,
} from 'lucide-react';
import confetti from 'canvas-confetti';

import CoverageChart from './components/CoverageChart';
import AmpliconChart from './components/AmpliconChart';
import GeneMap from './components/GeneMap';
import MutationMatrix from './components/MutationMatrix';
import VariantTable from './components/VariantTable';
import PhylogeneticTree from './components/PhylogeneticTree';
import ErrorBoundary from './components/ErrorBoundary';

import {
  computeStableReportId,
  formatConsequences,
  formatMetric,
  getBestReference,
  getCoDetections,
  getContigConsensus,
  getContigCoverage,
  getEvalBreadth,
  getEvalDepth,
  getEvalNRate,
  getMappedReads,
  getQcThresholds,
  getReferencesList,
  getRunConfig,
  getSampleStatus,
  getSampleStatusForRef,
  getSoftwareVersions,
  getWholeBreadth,
  validateUploadJSON,
} from './utils/reportSchema';
import {
  formatCohortMultiFasta,
  formatSampleFasta,
  triggerBrowserDownload,
} from './utils/exporters';

function StatusBadge({ status, reasons, basis, depth, wholeBreadth, evalBreadth, evalNRate, thresholds }) {
  const pillStyle =
    status === 'PASS'
      ? 'text-emerald-700 bg-emerald-50 border-emerald-200'
      : status === 'WARN'
      ? 'text-amber-700 bg-amber-50 border-amber-200'
      : 'text-red-700 bg-red-50 border-red-200';

  const IconComp = status === 'PASS' ? CheckCircle2 : status === 'WARN' ? AlertTriangle : XCircle;
  const iconColor =
    status === 'PASS' ? 'text-emerald-400' : status === 'WARN' ? 'text-amber-400' : 'text-red-400';

  return (
    <span className={`relative group cursor-help inline-flex items-center gap-1 text-[10px] font-bold px-2.5 py-0.5 rounded-full border select-none ${pillStyle}`}>
      <IconComp className="w-3 h-3" /> {status}
      <div
        className="absolute right-full top-1/2 -translate-y-1/2 mr-3 hidden group-hover:block z-50 w-64 p-3 bg-slate-900 border border-slate-700 text-slate-200 rounded-xl shadow-2xl text-left font-normal leading-relaxed select-text"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="font-bold text-white mb-2 text-xs tracking-wider flex items-center justify-between gap-2 border-b border-slate-800 pb-1.5">
          <span className="flex items-center gap-1.5">
            <IconComp className={`w-4 h-4 ${iconColor}`} />
            {status}
          </span>
          <span className="text-[9px] font-mono text-slate-400">{basis || 'whole_contig'}</span>
        </div>
        <div className="text-[10px] space-y-1">
          <p className="text-slate-300">
            Depth: <strong className="text-white">{depth.toLocaleString()}x</strong> · Callable Breadth:{' '}
            <strong className="text-white">{evalBreadth.toFixed(1)}%</strong>
            {basis !== 'whole_contig' && (
              <span className="text-slate-400"> (whole: {wholeBreadth.toFixed(1)}%)</span>
            )}{' '}
            · Internal N: <strong className="text-white">{evalNRate.toFixed(1)}%</strong>
          </p>
          {reasons && reasons.length > 0 && (
            <ul className="list-disc pl-4 space-y-0.5 text-slate-300 pt-1">
              {reasons.map((r, i) => (
                <li key={i} className="font-mono text-[9.5px]">
                  {r}
                </li>
              ))}
            </ul>
          )}
          <p className="text-[9px] text-slate-400 pt-1 border-t border-slate-800">
            PASS: &ge;{thresholds.pass_depth}x, &ge;{thresholds.pass_breadth}%, &le;{thresholds.pass_max_n}% N
          </p>
        </div>
      </div>
    </span>
  );
}

export default function App() {
  const [reportData, setReportData] = useState(null);
  const [selectedSampleId, setSelectedSampleId] = useState('');
  const [activeTab, setActiveTab] = useState('overview');
  const [isDragging, setIsDragging] = useState(false);
  const [uploadError, setUploadError] = useState('');
  const [sortField, setSortField] = useState('sample_id');
  const [sortAsc, setSortAsc] = useState(true);
  const [selectedRefName, setSelectedRefName] = useState('');

  useEffect(() => {
    fetch('/mock_run_report.json')
      .then((res) => {
        if (!res.ok) throw new Error('Demo report not found');
        return res.json();
      })
      .then((data) => {
        setReportData(data);
        if (data.samples && data.samples.length > 0) {
          setSelectedSampleId(data.samples[0].sample_id);
        }
        const refs = getReferencesList(data);
        if (refs.length > 0) {
          setSelectedRefName(refs[0].name);
        }
      })
      .catch((err) => {
        console.warn('Could not load demo report automatically:', err);
      });
  }, []);

  const celebrateSuccess = () => {
    confetti({
      particleCount: 80,
      spread: 60,
      origin: { y: 0.8 },
      colors: ['#DC2626', '#7C3AED', '#10B981'],
    });
  };

  const handleJSONData = (text) => {
    try {
      const data = JSON.parse(text);
      validateUploadJSON(data);
      setReportData(data);
      if (data.samples.length > 0) {
        const firstSample = data.samples[0];
        setSelectedSampleId(firstSample.sample_id);
        const bestRef = getBestReference(firstSample, data);
        const refs = getReferencesList(data);
        setSelectedRefName(bestRef ? bestRef.name : refs[0]?.name || '');
      }
      setUploadError('');
      setActiveTab('overview');
      celebrateSuccess();
    } catch (err) {
      setUploadError(`Failed to load JSON: ${err.message}`);
    }
  };

  const handleFileUpload = (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (event) => handleJSONData(event.target.result);
    reader.readAsText(file);
  };

  const handleDragOver = (e) => {
    e.preventDefault();
    setIsDragging(true);
  };

  const handleDragLeave = () => {
    setIsDragging(false);
  };

  const handleDrop = (e) => {
    e.preventDefault();
    setIsDragging(false);
    const file = e.dataTransfer.files[0];
    const isJsonFile =
      file &&
      (file.type === 'application/json' ||
        file.type === 'text/json' ||
        file.name?.toLowerCase().endsWith('.json'));
    if (isJsonFile) {
      const reader = new FileReader();
      reader.onload = (event) => handleJSONData(event.target.result);
      reader.readAsText(file);
    } else {
      setUploadError('Only .json run report files are supported.');
    }
  };

  if (!reportData) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen p-6 bg-slate-50">
        <div className="w-full max-w-lg p-8 glass-panel rounded-3xl flex flex-col gap-6 text-center animate-pulse">
          <Dna className="w-16 h-16 mx-auto text-red-600 animate-bounce" />
          <h2 className="text-2xl font-bold text-slate-900">Loading Viromics Suite...</h2>
          <p className="text-xs text-slate-500">Initializing analytical frameworks and visuals.</p>
        </div>
      </div>
    );
  }

  const minCovThreshold = reportData.min_coverage_threshold ?? 10;
  const runConfig = getRunConfig(reportData);
  const qcThresholds = getQcThresholds(reportData);
  const softwareVersions = getSoftwareVersions(reportData);
  const refsList = getReferencesList(reportData);
  const activeSample = reportData.samples.find((s) => s.sample_id === selectedSampleId) || reportData.samples[0];

  const handleSort = (field) => {
    if (sortField === field) {
      setSortAsc(!sortAsc);
    } else {
      setSortField(field);
      setSortAsc(true);
    }
  };

  const sortedSamples = [...reportData.samples].sort((a, b) => {
    if (!sortField) return 0;
    const refA = getBestReference(a, reportData)?.name;
    const refB = getBestReference(b, reportData)?.name;

    let valA;
    let valB;
    if (sortField === 'sample_id') {
      valA = a.sample_id;
      valB = b.sample_id;
    } else if (sortField === 'trimmed_reads') {
      valA = a.qc?.trimmed_reads ?? 0;
      valB = b.qc?.trimmed_reads ?? 0;
    } else if (sortField === 'mapped_reads') {
      valA = getMappedReads(a, refA) ?? -1;
      valB = getMappedReads(b, refB) ?? -1;
    } else if (sortField === 'mean_depth') {
      valA = getEvalDepth(a, refA);
      valB = getEvalDepth(b, refB);
    } else if (sortField === 'coverage_breadth') {
      valA = getEvalBreadth(a, refA);
      valB = getEvalBreadth(b, refB);
    } else if (sortField === 'n_content') {
      valA = getEvalNRate(a, refA);
      valB = getEvalNRate(b, refB);
    } else if (sortField === 'consensus_length') {
      valA = getContigConsensus(a, refA).length ?? 0;
      valB = getContigConsensus(b, refB).length ?? 0;
    } else if (sortField === 'snps_called') {
      valA = a.variants ? a.variants.filter((v) => !refA || !v.reference || v.reference === refA).length : 0;
      valB = b.variants ? b.variants.filter((v) => !refB || !v.reference || v.reference === refB).length : 0;
    } else if (sortField === 'status') {
      valA = getSampleStatus(a, reportData).status;
      valB = getSampleStatus(b, reportData).status;
    }

    if (typeof valA === 'string') {
      valA = valA.toLowerCase();
      valB = String(valB).toLowerCase();
    }
    if (valA < valB) return sortAsc ? -1 : 1;
    if (valA > valB) return sortAsc ? 1 : -1;
    return 0;
  });

  const handleDownloadCohortFasta = () => {
    const targetRef = selectedRefName || refsList[0]?.name || '';
    const fasta = formatCohortMultiFasta(reportData.samples, targetRef);
    if (!fasta) return;
    triggerBrowserDownload(fasta, `${reportData.run_id}_${targetRef || 'cohort'}_consensus.fasta`, 'text/plain;charset=utf-8;');
  };

  const passingSamplesCount = reportData.samples.filter((s) => getSampleStatus(s, reportData).status === 'PASS').length;
  const warnSamplesCount = reportData.samples.filter((s) => getSampleStatus(s, reportData).status === 'WARN').length;
  const failSamplesCount = reportData.samples.filter((s) => getSampleStatus(s, reportData).status === 'FAIL').length;
  const passPct = reportData.samples.length > 0 ? Math.round((passingSamplesCount / reportData.samples.length) * 100) : 0;

  const platformLabel = runConfig.platform ? String(runConfig.platform).toUpperCase() : 'UNKNOWN';
  const versionBadges = Object.entries(softwareVersions).slice(0, 4);

  return (
    <div className="min-h-screen pb-16 flex flex-col bg-slate-50">
      <header className="no-print app-header sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-6 h-16 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-9 h-9 rounded-xl bg-gradient-to-tr from-red-600 to-violet-600 flex items-center justify-center">
              <Dna className="w-5 h-5 text-white" />
            </div>
            <div>
              <h1 className="text-base font-bold text-white tracking-wide flex items-center gap-1.5">
                VIROMICS
                <span className="text-[10px] py-0.5 px-1.5 rounded bg-red-600/20 text-red-400 font-semibold font-mono">
                  v{reportData.schema_version || '1.0'}
                </span>
              </h1>
              <p className="text-[10px] text-slate-400">Genomic Epidemiology Suite</p>
            </div>
          </div>

          <div className="flex items-center gap-3">
            <div className="hidden md:flex items-center gap-2 bg-slate-800 border border-slate-700 px-3 py-1.5 rounded-xl text-xs">
              <Zap className="w-3.5 h-3.5 text-amber-400" />
              <span className="text-slate-400">Active Run:</span>
              <strong className="text-white font-mono">{reportData.run_id}</strong>
            </div>

            <label className="flex items-center gap-1.5 px-3 py-1.5 bg-slate-800 border border-slate-700 text-xs text-slate-100 font-semibold rounded-xl hover:bg-slate-700 cursor-pointer transition-all">
              <UploadCloud className="w-4 h-4 text-red-400" />
              Upload Run JSON
              <input type="file" accept=".json,application/json" onChange={handleFileUpload} className="hidden" />
            </label>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-6 w-full flex-1 flex flex-col gap-6 mt-8">
        {isDragging ? (
          <div
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleDrop}
            className="border-2 border-dashed border-red-500 bg-red-50 rounded-3xl p-16 flex flex-col items-center justify-center text-center gap-4"
          >
            <UploadCloud className="w-16 h-16 text-red-600 animate-bounce" />
            <h2 className="text-xl font-bold text-slate-900">Drop the report JSON file here</h2>
            <p className="text-sm text-slate-500">Instantly populate coverage, mutations, amplicons, and phylogeny.</p>
          </div>
        ) : (
          <div onDragOver={handleDragOver} className="flex flex-col gap-6">
            {uploadError && (
              <div className="bg-red-50 border border-red-200 text-red-700 text-xs p-3 rounded-xl flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 shrink-0" />
                <span>{uploadError}</span>
              </div>
            )}

            <div className="flex flex-wrap border-b border-slate-200 gap-1 sm:gap-2 no-print">
              {[
                { id: 'overview', label: 'Run Overview', icon: Grid },
                { id: 'details', label: 'Sample Details', icon: Sliders },
                { id: 'matrix', label: 'Mutational Sharing', icon: Layers },
                { id: 'tree', label: 'Phylogenetic Tree', icon: GitBranch },
                { id: 'report', label: 'Genomics Report', icon: FileText },
              ].map((tab) => {
                const TabIcon = tab.icon;
                const isActive = activeTab === tab.id;
                return (
                  <button
                    key={tab.id}
                    onClick={() => setActiveTab(tab.id)}
                    className={`flex items-center gap-2 px-4 py-2.5 text-xs font-semibold rounded-t-xl transition-all border-b-2 ${
                      isActive
                        ? 'border-red-600 text-red-600 bg-red-50/70'
                        : 'border-transparent text-slate-600 hover:text-slate-900'
                    }`}
                  >
                    <TabIcon className="w-4 h-4" />
                    {tab.label}
                  </button>
                );
              })}
            </div>

            {activeTab === 'overview' && (
              <ErrorBoundary componentName="Run Overview">
                <div className="flex flex-col gap-6">
                  <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
                    <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                      <div>
                        <span className="text-[10px] text-slate-500 font-bold uppercase tracking-wider">Total Samples</span>
                        <h3 className="text-2xl font-bold text-slate-900 mt-1">{reportData.samples.length}</h3>
                        <span className="text-[10px] text-slate-500 mt-1 block">Sequencing cohort size</span>
                      </div>
                      <div className="w-10 h-10 rounded-xl bg-red-50 flex items-center justify-center">
                        <BarChart3 className="w-5 h-5 text-red-600" />
                      </div>
                    </div>

                    <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                      <div>
                        <span className="text-[10px] text-slate-500 font-bold uppercase tracking-wider">Quality Pass Rate</span>
                        <h3 className="text-2xl font-bold text-emerald-600 mt-1">{passPct}%</h3>
                        <span className="text-[10px] text-slate-500 mt-1 block">
                          {passingSamplesCount} PASS, {warnSamplesCount} WARN, {failSamplesCount} FAIL
                        </span>
                      </div>
                      <div className="w-10 h-10 rounded-xl bg-emerald-50 flex items-center justify-center">
                        <ShieldCheck className="w-5 h-5 text-emerald-600" />
                      </div>
                    </div>

                    <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                      <div>
                        <span className="text-[10px] text-slate-500 font-bold uppercase tracking-wider">Reference Target</span>
                        <h3 className="text-sm font-bold text-slate-900 mt-1 font-mono truncate max-w-[165px]">
                          {refsList.length === 1
                            ? refsList[0].name
                            : refsList.length > 1
                            ? `${refsList.length} Reference Strains`
                            : 'Unknown'}
                        </h3>
                        <span className="text-[10px] text-slate-500 mt-1.5 block">
                          Length: {refsList.reduce((acc, r) => acc + (r.length || 0), 0).toLocaleString()} bp
                        </span>
                      </div>
                      <div className="w-10 h-10 rounded-xl bg-indigo-50 flex items-center justify-center">
                        <Dna className="w-5 h-5 text-indigo-600" />
                      </div>
                    </div>

                    <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                      <div className="min-w-0">
                        <span className="text-[10px] text-slate-500 font-bold uppercase tracking-wider">Run Configuration</span>
                        <h3 className="text-sm font-bold text-red-600 mt-1 font-mono">
                          {platformLabel} · &ge;{minCovThreshold}x · AF&ge;{runConfig.min_freq !== undefined ? `${(runConfig.min_freq * 100).toFixed(0)}%` : 'N/A'}
                        </h3>
                        <span className="text-[10px] text-slate-500 mt-1 block truncate">
                          {runConfig.use_iupac ? 'IUPAC Consensus' : `Haploid (AF>=${((runConfig.min_consensus_freq ?? 0.5) * 100).toFixed(0)}%)`}
                          {versionBadges.length > 0 ? ` · ${versionBadges.map(([k, v]) => `${k} ${v}`).join(', ')}` : ''}
                        </span>
                      </div>
                      <div className="w-10 h-10 rounded-xl bg-red-50 flex items-center justify-center shrink-0">
                        <Activity className="w-5 h-5 text-red-600" />
                      </div>
                    </div>
                  </div>

                  <div className="p-6 glass-panel rounded-2xl flex flex-col gap-4">
                    <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3">
                      <div>
                        <h3 className="text-lg font-semibold text-slate-900">Consensus Summary &amp; QC</h3>
                        <p className="text-xs text-slate-500">
                          Tier classifications evaluated dynamically from pipeline thresholds and callable-target intervals
                        </p>
                      </div>

                      {reportData.samples.length > 0 && (
                        <button
                          onClick={handleDownloadCohortFasta}
                          className="self-start sm:self-auto flex items-center gap-1.5 px-3 py-1.5 bg-white border border-slate-300 text-slate-700 text-xs font-semibold rounded-xl hover:bg-slate-50 transition-colors"
                        >
                          <Download className="w-3.5 h-3.5 text-red-600" />
                          Download Cohort FASTA
                        </button>
                      )}
                    </div>

                    <div className="flex flex-wrap gap-x-6 gap-y-1.5 text-[11px] text-slate-600 border-t border-slate-200 pt-2.5">
                      <span className="flex items-center gap-1.5">
                        <span className="w-2 h-2 rounded-full bg-emerald-500" />
                        <strong className="text-slate-900">PASS:</strong> &ge;{qcThresholds.pass_depth}x depth, &ge;{qcThresholds.pass_breadth}% breadth, &le;{qcThresholds.pass_max_n}% internal N
                      </span>
                      <span className="flex items-center gap-1.5">
                        <span className="w-2 h-2 rounded-full bg-amber-500" />
                        <strong className="text-slate-900">WARN:</strong> &ge;{qcThresholds.warn_depth}x depth, &ge;{qcThresholds.warn_breadth}% breadth, &le;{qcThresholds.warn_max_n}% internal N
                      </span>
                      <span className="flex items-center gap-1.5">
                        <span className="w-2 h-2 rounded-full bg-red-500" />
                        <strong className="text-slate-900">FAIL:</strong> &lt;{qcThresholds.warn_depth}x depth, &lt;{qcThresholds.warn_breadth}% breadth, or &gt;{qcThresholds.warn_max_n}% N
                      </span>
                    </div>

                    {reportData.samples.length === 0 ? (
                      <div className="p-8 text-center text-slate-500 border border-slate-200 rounded-xl bg-slate-50">
                        No samples present in this run report.
                      </div>
                    ) : (
                      <div className="overflow-x-auto border border-slate-200 rounded-xl bg-white max-w-full">
                        <table className="w-full text-left text-xs border-collapse">
                          <thead>
                            <tr className="border-b border-slate-200 text-slate-700 bg-slate-100 select-none">
                              {[
                                { id: 'sample_id', label: 'Sample ID' },
                                { id: 'trimmed_reads', label: 'Filtered Reads' },
                                { id: 'mapped_reads', label: 'Mapped Reads' },
                                { id: 'mean_depth', label: 'Mean Depth' },
                                { id: 'coverage_breadth', label: `Cov Breadth (≥${minCovThreshold}x)` },
                                { id: 'n_content', label: 'Internal N%' },
                                { id: 'consensus_length', label: 'Consensus Len' },
                                { id: 'snps_called', label: 'Variants' },
                              ].map((col) => (
                                <th
                                  key={col.id}
                                  className="p-3.5 font-semibold cursor-pointer hover:text-slate-900 transition-colors"
                                  onClick={() => handleSort(col.id)}
                                >
                                  <div className="flex items-center gap-1">
                                    {col.label} <ArrowUpDown className="w-3 h-3 text-red-600" />
                                  </div>
                                </th>
                              ))}
                              <th
                                className="p-3.5 font-semibold cursor-pointer hover:text-slate-900 transition-colors text-center"
                                onClick={() => handleSort('status')}
                              >
                                <div className="flex items-center justify-center gap-1">
                                  Status <ArrowUpDown className="w-3 h-3 text-red-600" />
                                </div>
                              </th>
                            </tr>
                          </thead>
                          <tbody>
                            {sortedSamples.map((sample) => {
                              const bestRef = getBestReference(sample, reportData);
                              const refName = bestRef?.name;
                              const cov = getContigCoverage(sample, refName);
                              const cons = getContigConsensus(sample, refName);
                              const statusInfo = getSampleStatus(sample, reportData);
                              const coDetections = getCoDetections(sample, reportData);
                              const mappedReads = getMappedReads(sample, refName);
                              const filteredReads = sample.qc?.trimmed_reads ?? 0;
                              const mappedPct =
                                mappedReads !== null && filteredReads > 0
                                  ? `${((mappedReads / filteredReads) * 100).toFixed(1)}%`
                                  : null;

                              const evalDepth = getEvalDepth(sample, refName);
                              const evalBreadth = getEvalBreadth(sample, refName);
                              const wholeBreadth = getWholeBreadth(sample, refName);
                              const evalNRate = getEvalNRate(sample, refName);
                              const primaryVarsCount = (sample.variants || []).filter(
                                (v) => !refName || !v.reference || v.reference === refName
                              ).length;

                              return (
                                <tr
                                  key={sample.sample_id}
                                  className="border-b border-slate-100 hover:bg-slate-50 cursor-pointer transition-colors"
                                  onClick={() => {
                                    setSelectedSampleId(sample.sample_id);
                                    if (refName) setSelectedRefName(refName);
                                    setActiveTab('details');
                                  }}
                                >
                                  <td className="p-3.5 font-semibold text-slate-900">
                                    <div>{sample.sample_id}</div>
                                    {coDetections.length > 1 ? (
                                      <div className="text-[10px] text-violet-700 font-bold mt-0.5">
                                        Co-detection: {coDetections.join(' + ')}
                                      </div>
                                    ) : (
                                      refName && <div className="text-[10px] text-slate-500 font-normal mt-0.5">{refName}</div>
                                    )}
                                  </td>
                                  <td className="p-3.5 font-mono text-slate-700">
                                    <div>{filteredReads.toLocaleString()}</div>
                                    <div className="text-[10px] text-slate-400">
                                      Dup: {formatMetric(sample.qc?.duplication_rate, { decimals: 1, suffix: '%', isRatioPercent: true })}
                                    </div>
                                  </td>
                                  <td className="p-3.5 font-mono text-indigo-700 font-medium">
                                    {mappedReads === null
                                      ? 'N/A'
                                      : `${mappedReads.toLocaleString()} (${mappedPct || '0.0%'})`}
                                  </td>
                                  <td className={`p-3.5 font-mono font-semibold ${
                                    evalDepth >= qcThresholds.pass_depth
                                      ? 'text-emerald-600'
                                      : evalDepth >= qcThresholds.warn_depth
                                      ? 'text-amber-600'
                                      : 'text-red-600'
                                  }`}>
                                    <div>{evalDepth.toLocaleString()}x</div>
                                    {cov.callable_mean_depth !== undefined && cov.mean_depth !== cov.callable_mean_depth && (
                                      <div className="text-[10px] text-slate-400 font-normal">
                                        whole: {cov.mean_depth.toLocaleString()}x
                                      </div>
                                    )}
                                  </td>
                                  <td className={`p-3.5 font-mono font-bold ${
                                    evalBreadth >= qcThresholds.pass_breadth
                                      ? 'text-emerald-600'
                                      : evalBreadth >= qcThresholds.warn_breadth
                                      ? 'text-amber-600'
                                      : 'text-red-600'
                                  }`}>
                                    <div>{evalBreadth.toFixed(2)}%</div>
                                    {cov.callable_breadth_min_cov !== undefined && cov.qc_basis !== 'whole_contig' && (
                                      <div className="text-[10px] text-slate-400 font-normal">
                                        whole: {wholeBreadth.toFixed(2)}%
                                      </div>
                                    )}
                                  </td>
                                  <td className={`p-3.5 font-mono font-bold ${
                                    evalNRate <= qcThresholds.pass_max_n
                                      ? 'text-emerald-600'
                                      : evalNRate <= qcThresholds.warn_max_n
                                      ? 'text-amber-600'
                                      : 'text-red-600'
                                  }`}>
                                    <div>{evalNRate.toFixed(2)}%</div>
                                    {cons.internal_n_content !== undefined && cov.qc_basis !== 'whole_contig' && (
                                      <div className="text-[10px] text-slate-400 font-normal">
                                        whole: {(cons.n_content ?? 0).toFixed(2)}%
                                      </div>
                                    )}
                                  </td>
                                  <td className="p-3.5 font-mono text-slate-700">
                                    {(cons.length ?? 0).toLocaleString()} bp
                                  </td>
                                  <td className="p-3.5 font-mono text-slate-700">{primaryVarsCount}</td>
                                  <td className="p-3.5 text-center">
                                    <StatusBadge
                                      status={statusInfo.status}
                                      reasons={statusInfo.reasons}
                                      basis={statusInfo.basis}
                                      depth={evalDepth}
                                      wholeBreadth={wholeBreadth}
                                      evalBreadth={evalBreadth}
                                      evalNRate={evalNRate}
                                      thresholds={qcThresholds}
                                    />
                                  </td>
                                </tr>
                              );
                            })}
                          </tbody>
                        </table>
                      </div>
                    )}
                  </div>
                </div>
              </ErrorBoundary>
            )}

            {activeTab === 'details' && activeSample && (
              <ErrorBoundary componentName="Sample Details">
                {(() => {
                  const targetRef = selectedRefName || getBestReference(activeSample, reportData)?.name || refsList[0]?.name || '';
                  const activeCov = getContigCoverage(activeSample, targetRef);
                  const activeCons = getContigConsensus(activeSample, targetRef);
                  const activeRefLength = refsList.find((r) => r.name === targetRef)?.length || 0;
                  const activeGenes = (targetRef && reportData.genes?.[targetRef]) || (Array.isArray(reportData.genes) ? reportData.genes : []);
                  const activeVariants = (activeSample.variants || []).filter(
                    (v) => !v.reference || v.reference === targetRef
                  );
                  const mappedReads = getMappedReads(activeSample, targetRef);
                  const filteredReads = activeSample.qc?.trimmed_reads || 0;
                  const mappedPct =
                    mappedReads !== null && filteredReads > 0
                      ? ((mappedReads / filteredReads) * 100).toFixed(1)
                      : null;
                  const evalDepth = getEvalDepth(activeSample, targetRef);
                  const evalBreadth = getEvalBreadth(activeSample, targetRef);

                  return (
                    <div className="flex flex-col gap-6">
                      <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-4 p-4 glass-panel rounded-2xl">
                        <div className="flex flex-wrap items-center gap-4">
                          <div className="flex items-center gap-2">
                            <span className="text-xs text-slate-500 font-bold uppercase tracking-wider">Active Sample:</span>
                            <select
                              value={selectedSampleId}
                              onChange={(e) => {
                                const newId = e.target.value;
                                setSelectedSampleId(newId);
                                const nextSample = reportData.samples.find((s) => s.sample_id === newId);
                                const bestRef = getBestReference(nextSample, reportData);
                                if (bestRef) setSelectedRefName(bestRef.name);
                              }}
                              className="bg-slate-50 border border-slate-200 text-slate-900 text-sm rounded-xl px-3.5 py-2 focus:outline-none cursor-pointer"
                            >
                              {[...reportData.samples]
                                .sort((a, b) => a.sample_id.localeCompare(b.sample_id, undefined, { numeric: true, sensitivity: 'base' }))
                                .map((s) => (
                                  <option key={s.sample_id} value={s.sample_id}>
                                    {s.sample_id}
                                  </option>
                                ))}
                            </select>
                          </div>

                          {refsList.length > 1 && (
                            <div className="flex items-center gap-2">
                              <span className="text-xs text-slate-500 font-bold uppercase tracking-wider">Reference:</span>
                              <select
                                value={targetRef}
                                onChange={(e) => setSelectedRefName(e.target.value)}
                                className="bg-slate-50 border border-slate-200 text-slate-900 text-sm rounded-xl px-3.5 py-2 focus:outline-none cursor-pointer"
                              >
                                {refsList.map((r) => (
                                  <option key={r.name} value={r.name}>
                                    {r.name} ({getEvalBreadth(activeSample, r.name).toFixed(1)}% cov)
                                  </option>
                                ))}
                              </select>
                            </div>
                          )}

                          <button
                            onClick={() => {
                              const fasta = formatSampleFasta(activeSample, targetRef);
                              if (!fasta) return;
                              triggerBrowserDownload(fasta, `${activeSample.sample_id}_${targetRef || 'consensus'}.fasta`, 'text/plain;charset=utf-8;');
                            }}
                            disabled={!activeCons?.sequence}
                            className="flex items-center gap-1.5 px-3 py-2 bg-white border border-slate-300 text-slate-700 text-xs font-semibold rounded-xl hover:bg-slate-50 disabled:opacity-40"
                          >
                            <Download className="w-3.5 h-3.5 text-red-600" />
                            Sample FASTA
                          </button>
                        </div>

                        <div className="flex flex-wrap items-center gap-x-3 gap-y-1 text-xs text-slate-600">
                          <span>
                            Mapped:{' '}
                            <strong className="text-indigo-700 font-mono">
                              {mappedReads === null
                                ? 'N/A'
                                : `${mappedReads.toLocaleString()} (${mappedPct || '0.0'}%)`}
                            </strong>
                          </span>
                          <span className="text-slate-300">|</span>
                          <span>
                            Callable Depth:{' '}
                            <strong className={`font-mono ${evalDepth >= qcThresholds.pass_depth ? 'text-emerald-600' : evalDepth >= qcThresholds.warn_depth ? 'text-amber-600' : 'text-red-600'}`}>
                              {evalDepth}x
                            </strong>
                          </span>
                          <span className="text-slate-300">|</span>
                          <span>
                            Breadth (&ge;{minCovThreshold}x):{' '}
                            <strong className="font-mono text-slate-900">{evalBreadth.toFixed(2)}%</strong>
                          </span>
                          <span className="text-slate-300">|</span>
                          <span>
                            Q30 Trim:{' '}
                            <strong className="font-mono text-slate-900">
                              {formatMetric(activeSample.qc?.q30_rate_trimmed, { decimals: 1, suffix: '%', isRatioPercent: true })}
                            </strong>
                          </span>
                          <span className="text-slate-300">|</span>
                          <span>
                            Dup:{' '}
                            <strong className="font-mono text-slate-900">
                              {formatMetric(activeSample.qc?.duplication_rate, { decimals: 1, suffix: '%', isRatioPercent: true })}
                            </strong>
                          </span>
                        </div>
                      </div>

                      {activeSample.metadata && Object.keys(activeSample.metadata).length > 0 && (
                        <div className="flex flex-wrap items-center gap-x-6 gap-y-2 text-xs bg-slate-100 border border-slate-200 rounded-2xl p-4">
                          {Object.entries(activeSample.metadata)
                            .filter(([k]) => k !== 'sample_id')
                            .map(([key, val]) => (
                              <span key={key}>
                                <span className="text-[10px] text-slate-500 uppercase font-bold tracking-wider mr-2">
                                  {key.replace(/_/g, ' ')}:
                                </span>
                                <strong className="text-slate-900 font-semibold">{val}</strong>
                              </span>
                            ))}
                        </div>
                      )}

                      <CoverageChart
                        coverageData={activeCov.coverage_profile}
                        minCoverageThreshold={minCovThreshold}
                        lowCoverageIntervals={activeCov.low_coverage_intervals}
                        callableTargetIntervals={activeCov.callable_target_intervals}
                        sampleId={activeSample.sample_id}
                        selectedRefName={targetRef}
                      />

                      <AmpliconChart
                        amplicons={activeSample.amplicons || []}
                        selectedRefName={targetRef}
                        sampleId={activeSample.sample_id}
                        minCoverageThreshold={minCovThreshold}
                      />

                      <GeneMap
                        genes={activeGenes}
                        variants={activeVariants}
                        refLength={activeRefLength}
                        sampleId={activeSample.sample_id}
                        selectedRefName={targetRef}
                      />

                      <VariantTable variants={activeVariants} sampleId={activeSample.sample_id} />
                    </div>
                  );
                })()}
              </ErrorBoundary>
            )}

            {activeTab === 'matrix' && (
              <ErrorBoundary componentName="Mutational Sharing">
                <MutationMatrix
                  samples={reportData.samples}
                  selectedRefName={selectedRefName}
                  setSelectedRefName={setSelectedRefName}
                  refs={refsList}
                  minCoverageThreshold={minCovThreshold}
                  runConfig={runConfig}
                />
              </ErrorBoundary>
            )}

            {activeTab === 'tree' && (
              <ErrorBoundary componentName="Phylogenetic Tree">
                <PhylogeneticTree
                  samples={reportData.samples}
                  selectedRefName={selectedRefName}
                  setSelectedRefName={setSelectedRefName}
                  refs={refsList}
                  reportData={reportData}
                  minCoverageThreshold={minCovThreshold}
                />
              </ErrorBoundary>
            )}

            {activeTab === 'report' && activeSample && (
              <ErrorBoundary componentName="Genomics Report">
                {(() => {
                  const bestRef = getBestReference(activeSample, reportData);
                  const pdfRefName = bestRef ? bestRef.name : refsList[0]?.name || '';
                  const pdfCov = getContigCoverage(activeSample, pdfRefName);
                  const pdfCons = getContigConsensus(activeSample, pdfRefName);
                  const coDetections = getCoDetections(activeSample, reportData);
                  const stableReportId = computeStableReportId(reportData.run_id, activeSample.sample_id);

                  return (
                    <div className="flex flex-col gap-6">
                      <div className="flex flex-wrap items-center justify-between gap-4 p-4 glass-panel rounded-2xl no-print">
                        <div>
                          <h4 className="text-sm font-semibold text-slate-900">Genomic Epidemiological PDF Generator</h4>
                          <p className="text-xs text-slate-500">Multi-page print layout with run configuration and software provenance.</p>
                        </div>

                        <div className="flex items-center gap-3">
                          <select
                            value={selectedSampleId}
                            onChange={(e) => setSelectedSampleId(e.target.value)}
                            className="bg-slate-50 border border-slate-200 text-slate-900 text-xs rounded-xl px-3 py-1.5 focus:outline-none cursor-pointer"
                          >
                            {[...reportData.samples]
                              .sort((a, b) => a.sample_id.localeCompare(b.sample_id, undefined, { numeric: true, sensitivity: 'base' }))
                              .map((s) => (
                                <option key={s.sample_id} value={s.sample_id}>
                                  {s.sample_id}
                                </option>
                              ))}
                          </select>

                          <button
                            onClick={() => window.print()}
                            className="flex items-center gap-1.5 px-3 py-1.5 bg-red-600 text-white text-xs font-bold rounded-xl hover:bg-red-700 active:scale-95 transition-all"
                          >
                            <Printer className="w-4 h-4" />
                            Print / Save PDF
                          </button>
                        </div>
                      </div>

                      <div className="p-8 bg-white text-slate-800 rounded-3xl max-w-4xl mx-auto shadow-lg flex flex-col gap-8 border border-slate-200 print-card">
                        <div className="flex items-start justify-between border-b border-slate-200 pb-6 print-avoid-break">
                          <div className="flex flex-col gap-1.5">
                            <span className="text-[10px] font-bold uppercase tracking-wider text-red-600">Genomic Epidemiology Report</span>
                            <h2 className="text-2xl font-bold text-slate-900 leading-none">Pathogen Genomics Summary</h2>
                            <span className="text-xs text-slate-500 font-mono mt-1">Report ID: {stableReportId}</span>
                          </div>
                          <div className="flex flex-col items-end text-right gap-1 text-xs">
                            <span className="font-semibold text-slate-700">Viromics ({platformLabel})</span>
                            <span className="text-slate-500">Run ID: {reportData.run_id}</span>
                            <span className="text-slate-400 font-mono text-[10px]">
                              Schema v{reportData.schema_version || '1.0'} · Depth&ge;{minCovThreshold}x
                            </span>
                          </div>
                        </div>

                        <div className="grid grid-cols-2 md:grid-cols-3 gap-6 bg-slate-50 p-5 rounded-2xl border border-slate-200 print-avoid-break">
                          <div className="flex flex-col gap-1">
                            <span className="text-[10px] uppercase font-bold text-slate-500">Sample Identifier</span>
                            <strong className="text-sm font-semibold text-slate-900">{activeSample.sample_id}</strong>
                          </div>
                          <div className="flex flex-col gap-1">
                            <span className="text-[10px] uppercase font-bold text-slate-500">Reference Target</span>
                            <strong className="text-sm font-semibold text-slate-900 font-mono">
                              {coDetections.length > 1 ? `Co-detection: ${coDetections.join(' + ')}` : pdfRefName}
                            </strong>
                          </div>
                          <div className="flex flex-col gap-1">
                            <span className="text-[10px] uppercase font-bold text-slate-500">Consensus Length</span>
                            <strong className="text-sm font-semibold text-slate-900 font-mono">
                              {coDetections.length > 1
                                ? coDetections
                                    .map((c) => `${getContigConsensus(activeSample, c).length ?? 0} bp (${c})`)
                                    .join(' | ')
                                : `${(pdfCons.length ?? 0).toLocaleString()} bp`}
                            </strong>
                          </div>
                        </div>

                        <div className="print-avoid-break">
                          <h3 className="text-xs font-bold text-slate-900 border-b border-slate-200 pb-2 mb-4 uppercase tracking-wider">
                            Sequencing Quality Assurance
                          </h3>
                          {coDetections.length > 1 ? (
                            <div className="overflow-hidden border border-slate-200 rounded-xl">
                              <table className="w-full text-left text-xs border-collapse">
                                <thead>
                                  <tr className="bg-slate-100 text-slate-700 font-semibold border-b border-slate-200">
                                    <th className="p-3">Detected Strain</th>
                                    <th className="p-3">Mapped Reads</th>
                                    <th className="p-3">Callable Depth</th>
                                    <th className="p-3">Callable Breadth (&ge;{minCovThreshold}x)</th>
                                    <th className="p-3">Internal N Rate</th>
                                    <th className="p-3 text-center">Status</th>
                                  </tr>
                                </thead>
                                <tbody>
                                  {coDetections.map((c) => {
                                    const statusInfo = getSampleStatusForRef(activeSample, c, reportData);
                                    const mr = getMappedReads(activeSample, c);
                                    const filtered = activeSample.qc?.trimmed_reads || 0;
                                    const pct = mr !== null && filtered > 0 ? ((mr / filtered) * 100).toFixed(1) : null;
                                    const depth = getEvalDepth(activeSample, c);
                                    const breadth = getEvalBreadth(activeSample, c);
                                    const nRate = getEvalNRate(activeSample, c);
                                    return (
                                      <tr key={c} className="border-b border-slate-100 text-slate-700">
                                        <td className="p-3 font-semibold text-slate-900">{c}</td>
                                        <td className="p-3 font-mono">{mr === null ? 'N/A' : `${mr.toLocaleString()} (${pct || '0.0'}%)`}</td>
                                        <td className="p-3 font-mono font-semibold">{depth.toLocaleString()}x</td>
                                        <td className="p-3 font-mono font-bold">{breadth.toFixed(2)}%</td>
                                        <td className="p-3 font-mono font-bold">{nRate.toFixed(2)}%</td>
                                        <td className="p-3 text-center font-bold">{statusInfo.status}</td>
                                      </tr>
                                    );
                                  })}
                                </tbody>
                              </table>
                            </div>
                          ) : (
                            <div className="grid grid-cols-2 md:grid-cols-5 gap-4 text-center">
                              <div className="p-3 border border-slate-200 rounded-xl bg-slate-50">
                                <span className="text-[9px] uppercase font-bold text-slate-500 block">Filtered Reads</span>
                                <strong className="text-sm text-slate-900 font-semibold font-mono">
                                  {(activeSample.qc?.trimmed_reads ?? 0).toLocaleString()}
                                </strong>
                              </div>
                              <div className="p-3 border border-slate-200 rounded-xl bg-slate-50">
                                <span className="text-[9px] uppercase font-bold text-slate-500 block">Mapped Reads</span>
                                <strong className="text-sm text-indigo-700 font-semibold font-mono">
                                  {(() => {
                                    const mr = getMappedReads(activeSample, pdfRefName);
                                    return mr === null ? 'N/A' : mr.toLocaleString();
                                  })()}
                                </strong>
                              </div>
                              <div className="p-3 border border-slate-200 rounded-xl bg-slate-50">
                                <span className="text-[9px] uppercase font-bold text-slate-500 block">Callable Depth</span>
                                <strong className="text-sm text-slate-900 font-semibold font-mono">
                                  {getEvalDepth(activeSample, pdfRefName).toLocaleString()}x
                                </strong>
                              </div>
                              <div className="p-3 border border-slate-200 rounded-xl bg-slate-50">
                                <span className="text-[9px] uppercase font-bold text-slate-500 block">Breadth (&ge;{minCovThreshold}x)</span>
                                <strong className="text-sm text-slate-900 font-bold font-mono">
                                  {getEvalBreadth(activeSample, pdfRefName).toFixed(2)}%
                                </strong>
                              </div>
                              <div className="p-3 border border-slate-200 rounded-xl bg-slate-50">
                                <span className="text-[9px] uppercase font-bold text-slate-500 block">Internal N Rate</span>
                                <strong className="text-sm text-slate-900 font-bold font-mono">
                                  {getEvalNRate(activeSample, pdfRefName).toFixed(2)}%
                                </strong>
                                {pdfCov.qc_basis && pdfCov.qc_basis !== 'whole_contig' && (
                                  <span className="text-[9px] text-slate-400 block">
                                    whole N: {(pdfCons.n_content ?? 0).toFixed(2)}%
                                  </span>
                                )}
                              </div>
                            </div>
                          )}
                        </div>

                        {Object.keys(softwareVersions).length > 0 && (
                          <div className="bg-slate-50 border border-slate-200 rounded-xl p-3 text-[10px] font-mono text-slate-600 print-avoid-break">
                            <strong className="text-slate-800 mr-2">Provenance:</strong>
                            {Object.entries(softwareVersions)
                              .map(([tool, ver]) => `${tool}=${ver}`)
                              .join(' · ')}
                          </div>
                        )}

                        <div className="flex flex-col gap-6">
                          {(() => {
                            const allVariants = activeSample.variants || [];
                            const renderMutationsTable = (titleLabel, filteredVars) => (
                              <div className="flex flex-col gap-3">
                                <div className="flex justify-between items-center border-b border-slate-200 pb-2">
                                  <h3 className="text-xs font-bold text-slate-900 uppercase tracking-wider">
                                    Identified Variants: {titleLabel}
                                  </h3>
                                  <span className="text-[10px] text-slate-600 font-semibold bg-slate-100 px-2 py-0.5 rounded-full">
                                    {filteredVars.length} variants
                                  </span>
                                </div>
                                <div className="overflow-hidden border border-slate-200 rounded-xl">
                                  <table className="w-full text-left text-xs border-collapse">
                                    <thead>
                                      <tr className="bg-slate-100 border-b border-slate-200 text-slate-600 font-semibold">
                                        <th className="p-3">Position</th>
                                        <th className="p-3">Allele</th>
                                        <th className="p-3">Gene Region</th>
                                        <th className="p-3">Protein Change</th>
                                        <th className="p-3">Type</th>
                                        <th className="p-3 text-right">Allele Freq</th>
                                      </tr>
                                    </thead>
                                    <tbody>
                                      {filteredVars.length > 0 ? (
                                        filteredVars.map((v, i) => {
                                          const cons = formatConsequences(v);
                                          return (
                                            <tr key={i} className="border-b border-slate-100 text-slate-700">
                                              <td className="p-3 font-mono font-medium">{v.position}</td>
                                              <td className="p-3 font-mono font-bold text-slate-900">
                                                {v.ref} &gt; {v.alt}
                                              </td>
                                              <td className="p-3">
                                                <span className="px-2 py-0.5 rounded text-[10px] font-semibold bg-indigo-50 text-indigo-700 border border-indigo-100">
                                                  {v.gene}
                                                </span>
                                              </td>
                                              <td className="p-3 font-mono text-[11px]">
                                                {cons.length > 0 ? (
                                                  cons.map((c, cIdx) => (
                                                    <span key={cIdx} className="mr-1.5">
                                                      {cons.length > 1 ? `${c.gene}: ` : ''}
                                                      <span className={c.kind === 'missense' || c.kind === 'nonsense' ? 'font-bold text-red-600' : c.kind === 'frameshift' || c.kind === 'indel' ? 'text-amber-700 font-semibold' : 'text-slate-500'}>
                                                        {c.label}
                                                      </span>
                                                    </span>
                                                  ))
                                                ) : (
                                                  <span className="text-slate-400">-</span>
                                                )}
                                              </td>
                                              <td className="p-3">{v.type}</td>
                                              <td className="p-3 font-mono text-right text-red-600 font-bold">
                                                {(v.frequency * 100).toFixed(1)}%
                                              </td>
                                            </tr>
                                          );
                                        })
                                      ) : (
                                        <tr>
                                          <td colSpan={6} className="p-4 text-center text-slate-400">
                                            No variants identified for this reference strain.
                                          </td>
                                        </tr>
                                      )}
                                    </tbody>
                                  </table>
                                </div>
                              </div>
                            );

                            if (coDetections.length > 1) {
                              return coDetections.map((c) => (
                                <div key={c} className="mt-2">
                                  {renderMutationsTable(
                                    c,
                                    allVariants.filter((v) => v.reference === c)
                                  )}
                                </div>
                              ));
                            }
                            return renderMutationsTable(
                              pdfRefName,
                              allVariants.filter((v) => !v.reference || v.reference === pdfRefName)
                            );
                          })()}
                        </div>

                        <div className="flex justify-between items-center border-t border-slate-200 pt-8 mt-4 text-xs print-avoid-break">
                          <div className="flex flex-col gap-1">
                            <div className="w-32 border-b border-slate-300 h-8" />
                            <span className="text-slate-400 uppercase font-bold text-[9px] tracking-wider">Clinician / Assayer Signature</span>
                          </div>
                          <div className="flex flex-col gap-1 text-right">
                            <div className="w-32 border-b border-slate-300 h-8 self-end" />
                            <span className="text-slate-400 uppercase font-bold text-[9px] tracking-wider">Lab Director Approval</span>
                          </div>
                        </div>
                      </div>
                    </div>
                  );
                })()}
              </ErrorBoundary>
            )}
          </div>
        )}
      </main>
    </div>
  );
}
