import React, { useState, useEffect } from 'react';
import {
  UploadCloud,
  Grid,
  Sliders,
  GitBranch,
  AlertTriangle,
  CheckCircle2,
  XCircle,
  FileText,
  Layers,
  Printer,
  ArrowUpDown,
  Download,
} from 'lucide-react';

import CoverageChart from './components/CoverageChart';
import AmpliconChart from './components/AmpliconChart';
import GeneMap from './components/GeneMap';
import MutationMatrix from './components/MutationMatrix';
import VariantTable from './components/VariantTable';
import PhylogeneticTree from './components/PhylogeneticTree';
import ErrorBoundary from './components/ErrorBoundary';

import {
  computeStableReportId,
  deduplicateConsequences,
  formatConsequences,
  formatDepth,
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

function ViromicsLogoMark() {
  return (
    <div className="w-8 h-8 rounded-[4px] bg-[#EAF3F5] border border-[#B9D6DD] flex items-center justify-center shrink-0">
      <svg viewBox="0 0 24 24" className="w-4 h-4" fill="none" aria-hidden="true">
        <path
          d="M5 7L12 18L19 7"
          stroke="#0F5C6E"
          strokeWidth="2.2"
          strokeLinecap="round"
          strokeLinejoin="round"
        />
        <circle cx="12" cy="6.5" r="1.75" fill="#0F5C6E" />
      </svg>
    </div>
  );
}

function getQcBadgeClasses(status) {
  if (status === 'PASS') return 'text-[#16654A] bg-[#ECF7F2] border-[#B5DFCE]';
  if (status === 'WARN') return 'text-[#8C5308] bg-[#FDF6E7] border-[#EBD29B]';
  return 'text-[#A62C22] bg-[#FCEEEB] border-[#EBBDB8]';
}

function getQcIconColor(status) {
  if (status === 'PASS') return 'text-[#16654A]';
  if (status === 'WARN') return 'text-[#8C5308]';
  return 'text-[#A62C22]';
}

function StatusBadge({ status, reasons, basis, depth, wholeBreadth, evalBreadth, evalNRate, thresholds }) {
  const pillStyle = getQcBadgeClasses(status);
  const IconComp = status === 'PASS' ? CheckCircle2 : status === 'WARN' ? AlertTriangle : XCircle;
  const iconColor = getQcIconColor(status);

  return (
    <span className={`relative group cursor-help inline-flex items-center gap-1 text-[11px] font-semibold px-2 py-0.5 rounded-[4px] border select-none tabular-nums ${pillStyle}`}>
      <IconComp className="w-3 h-3" /> {status}
      <div
        className="absolute right-full top-1/2 -translate-y-1/2 mr-2 hidden group-hover:block z-50 w-64 p-3 bg-white border border-[#CBD2D9] text-[#111827] rounded-[6px] shadow-tooltip text-left font-normal leading-relaxed select-text"
        onClick={(e) => e.stopPropagation()}
      >
        <div className="font-semibold text-[#111827] mb-2 text-xs flex items-center justify-between gap-2 border-b border-[#E2E6EA] pb-1.5">
          <span className="flex items-center gap-1.5">
            <IconComp className={`w-3.5 h-3.5 ${iconColor}`} />
            {status}
          </span>
          <span className="text-[10px] font-mono text-[#6B7280]">{basis || 'whole_contig'}</span>
        </div>
        <div className="text-[11px] space-y-1 text-[#4B5563]">
          <p>
            Depth: <strong className="text-[#111827] font-mono">{formatDepth(depth)}</strong> · Callable Breadth:{' '}
            <strong className="text-[#111827] font-mono">{evalBreadth.toFixed(1)}%</strong>
            {basis !== 'whole_contig' && (
              <span className="text-[#6B7280]"> (whole: {wholeBreadth.toFixed(1)}%)</span>
            )}{' '}
            · Internal N: <strong className="text-[#111827] font-mono">{evalNRate.toFixed(1)}%</strong>
          </p>
          {reasons && reasons.length > 0 && (
            <ul className="list-disc pl-4 space-y-0.5 text-[#4B5563] pt-1">
              {reasons.map((r, i) => (
                <li key={i} className="font-mono text-[10px]">
                  {r}
                </li>
              ))}
            </ul>
          )}
          <p className="text-[10px] text-[#6B7280] pt-1.5 border-t border-[#E2E6EA] font-mono">
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
      <div className="flex flex-col items-center justify-center min-h-screen p-6 bg-[#F8F9FA]">
        <div className="w-full max-w-md p-8 glass-panel flex flex-col items-center gap-3 text-center">
          <ViromicsLogoMark />
          <h2 className="text-[18px] font-semibold text-[#111827]">Loading Viromics Suite</h2>
          <p className="text-[13px] text-[#6B7280]">Initializing cohort coverage, variants, and phylogenetics.</p>
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
  const versionEntries = Object.entries(softwareVersions);

  return (
    <div className="min-h-screen pb-14 flex flex-col bg-[#F8F9FA]">
      <header className="no-print app-header sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-6 h-14 flex items-center justify-between gap-4">
          <div className="flex items-center gap-3">
            <ViromicsLogoMark />
            <div className="flex flex-col">
              <div className="flex items-center gap-2">
                <h1 className="text-[15px] font-bold text-[#111827] tracking-tight leading-none">
                  VIROMICS
                </h1>
                <span className="text-[10px] py-0.5 px-1.5 rounded-[4px] bg-[#EAF3F5] text-[#0B4654] border border-[#B9D6DD] font-semibold font-mono leading-none">
                  v{reportData.schema_version || '1.0'}
                </span>
              </div>
              <p className="text-[11px] text-[#6B7280] leading-tight mt-0.5">
                Genomic Epidemiology Suite
              </p>
            </div>
          </div>

          <div className="flex items-center gap-2.5">
            <div className="hidden md:flex items-center gap-2 bg-[#F8F9FA] border border-[#E2E6EA] px-2.5 py-1.5 rounded-[4px] text-xs">
              <span className="text-[#6B7280]">Run:</span>
              <strong className="text-[#111827] font-mono font-semibold">{reportData.run_id}</strong>
            </div>

            <label className="flex items-center gap-1.5 px-3 py-1.5 bg-[#0F5C6E] hover:bg-[#0B4654] text-xs text-white font-semibold rounded-[4px] cursor-pointer transition-colors">
              <UploadCloud className="w-3.5 h-3.5" />
              Upload Run JSON
              <input type="file" accept=".json,application/json" onChange={handleFileUpload} className="hidden" />
            </label>
          </div>
        </div>
      </header>

      <main className="max-w-7xl mx-auto px-6 w-full flex-1 flex flex-col gap-5 mt-6">
        {isDragging ? (
          <div
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleDrop}
            className="border-2 border-dashed border-[#0F5C6E] bg-[#EAF3F5] rounded-[6px] p-16 flex flex-col items-center justify-center text-center gap-3"
          >
            <UploadCloud className="w-10 h-10 text-[#0F5C6E]" />
            <h2 className="text-[18px] font-semibold text-[#111827]">Drop the report JSON file here</h2>
            <p className="text-[13px] text-[#4B5563]">Populate cohort coverage, mutations, amplicons, and phylogeny.</p>
          </div>
        ) : (
          <div onDragOver={handleDragOver} className="flex flex-col gap-5">
            {uploadError && (
              <div className="bg-[#FCEEEB] border border-[#EBBDB8] text-[#A62C22] text-xs px-3.5 py-2.5 rounded-[4px] flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 shrink-0" />
                <span>{uploadError}</span>
              </div>
            )}

            <nav className="flex flex-wrap border-b border-[#E2E6EA] gap-1 no-print" aria-label="Dashboard sections">
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
                    className={`flex items-center gap-1.5 px-3.5 py-2.5 text-xs font-semibold transition-colors border-b-2 -mb-px ${
                      isActive
                        ? 'border-[#0F5C6E] text-[#0F5C6E]'
                        : 'border-transparent text-[#4B5563] hover:text-[#111827]'
                    }`}
                  >
                    <TabIcon className="w-3.5 h-3.5" />
                    {tab.label}
                  </button>
                );
              })}
            </nav>

            {activeTab === 'overview' && (
              <ErrorBoundary componentName="Run Overview">
                <div className="flex flex-col gap-4">
                  <div className="glass-panel grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 divide-y sm:divide-y-0 sm:divide-x divide-[#E2E6EA]">
                    <div className="px-5 py-4 flex flex-col justify-between">
                      <span className="text-[11px] text-[#6B7280] font-semibold uppercase tracking-[0.06em]">
                        Total Samples
                      </span>
                      <div className="text-[22px] font-semibold text-[#111827] font-mono tabular-nums mt-1 leading-tight">
                        {reportData.samples.length}
                      </div>
                      <span className="text-[12px] text-[#6B7280] mt-1">Sequencing cohort size</span>
                    </div>

                    <div className="px-5 py-4 flex flex-col justify-between">
                      <span className="text-[11px] text-[#6B7280] font-semibold uppercase tracking-[0.06em]">
                        Quality Pass Rate
                      </span>
                      <div
                        className={`text-[22px] font-semibold font-mono tabular-nums mt-1 leading-tight ${
                          passPct === 100
                            ? 'text-[#16654A]'
                            : failSamplesCount > 0
                            ? 'text-[#111827]'
                            : 'text-[#8C5308]'
                        }`}
                      >
                        {passPct}%
                      </div>
                      <span className="text-[12px] text-[#4B5563] mt-1 tabular-nums">
                        <strong className="text-[#16654A] font-semibold">{passingSamplesCount} PASS</strong>
                        {' · '}
                        <strong className="text-[#8C5308] font-semibold">{warnSamplesCount} WARN</strong>
                        {' · '}
                        <strong className="text-[#A62C22] font-semibold">{failSamplesCount} FAIL</strong>
                      </span>
                    </div>

                    <div className="px-5 py-4 flex flex-col justify-between min-w-0">
                      <span className="text-[11px] text-[#6B7280] font-semibold uppercase tracking-[0.06em]">
                        Reference Target
                      </span>
                      <div
                        className="text-[13px] font-semibold text-[#111827] mt-1.5 font-mono truncate"
                        title={refsList.map((r) => r.name).join(', ')}
                      >
                        {refsList.length === 1
                          ? refsList[0].name
                          : refsList.length > 1
                          ? `${refsList.length} Reference Strains`
                          : 'Unknown'}
                      </div>
                      <span className="text-[12px] text-[#6B7280] mt-1 font-mono tabular-nums">
                        Length: {refsList.reduce((acc, r) => acc + (r.length || 0), 0).toLocaleString()} bp
                      </span>
                    </div>

                    <div className="px-5 py-4 flex flex-col justify-between min-w-0">
                      <span className="text-[11px] text-[#6B7280] font-semibold uppercase tracking-[0.06em]">
                        Run Configuration
                      </span>
                      <div className="text-[13px] font-semibold text-[#0F5C6E] mt-1.5 font-mono tabular-nums">
                        {platformLabel} · &ge;{minCovThreshold}x · AF&ge;
                        {runConfig.min_freq !== undefined ? `${(runConfig.min_freq * 100).toFixed(0)}%` : 'N/A'}
                      </div>
                      <span className="text-[12px] text-[#6B7280] mt-1 font-mono tabular-nums">
                        {runConfig.use_iupac
                          ? 'IUPAC Consensus'
                          : `Haploid (AF>=${((runConfig.min_consensus_freq ?? 0.5) * 100).toFixed(0)}%)`}
                      </span>
                    </div>
                  </div>

                  {versionEntries.length > 0 && (
                    <div className="px-4 py-2 bg-[#F1F3F5] border border-[#E2E6EA] rounded-[6px] flex flex-wrap items-center gap-x-4 gap-y-1 text-[11px] font-mono text-[#4B5563]">
                      <span className="font-sans font-semibold uppercase tracking-[0.06em] text-[10px] text-[#6B7280]">
                        Pipeline Versions
                      </span>
                      {versionEntries.map(([tool, ver]) => (
                        <span key={tool} className="whitespace-nowrap">
                          <strong className="text-[#111827] font-semibold">{tool}</strong> {ver}
                        </span>
                      ))}
                    </div>
                  )}

                  <div className="p-5 glass-panel flex flex-col gap-4">
                    <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-3">
                      <div>
                        <h3 className="text-[16px] font-semibold text-[#111827]">Consensus Summary &amp; QC</h3>
                        <p className="text-[12px] text-[#6B7280] mt-0.5">
                          Tier classifications evaluated dynamically from pipeline thresholds and callable-target intervals
                        </p>
                      </div>

                      {reportData.samples.length > 0 && (
                        <button
                          onClick={handleDownloadCohortFasta}
                          className="self-start sm:self-auto flex items-center gap-1.5 px-3 py-1.5 bg-white border border-[#CBD2D9] text-[#111827] text-xs font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] transition-colors whitespace-nowrap"
                        >
                          <Download className="w-3.5 h-3.5 text-[#0F5C6E]" />
                          Download Cohort FASTA
                        </button>
                      )}
                    </div>

                    <div className="flex flex-wrap gap-x-6 gap-y-1.5 text-[12px] text-[#4B5563] border-t border-[#E2E6EA] pt-2.5 tabular-nums">
                      <span className="flex items-center gap-1.5">
                        <span className="w-2 h-2 rounded-full bg-[#2A8C6A]" />
                        <strong className="text-[#16654A]">PASS:</strong> &ge;{qcThresholds.pass_depth}x depth, &ge;{qcThresholds.pass_breadth}% breadth, &le;{qcThresholds.pass_max_n}% internal N
                      </span>
                      <span className="flex items-center gap-1.5">
                        <span className="w-2 h-2 rounded-full bg-[#CC8400]" />
                        <strong className="text-[#8C5308]">WARN:</strong> &ge;{qcThresholds.warn_depth}x depth, &ge;{qcThresholds.warn_breadth}% breadth, &le;{qcThresholds.warn_max_n}% internal N
                      </span>
                      <span className="flex items-center gap-1.5">
                        <span className="w-2 h-2 rounded-full bg-[#C0392B]" />
                        <strong className="text-[#A62C22]">FAIL:</strong> &lt;{qcThresholds.warn_depth}x depth, &lt;{qcThresholds.warn_breadth}% breadth, or &gt;{qcThresholds.warn_max_n}% N
                      </span>
                    </div>

                    {reportData.samples.length === 0 ? (
                      <div className="p-8 text-center text-[#6B7280] border border-[#E2E6EA] rounded-[6px] bg-[#F8F9FA] text-[13px]">
                        No samples present in this run report.
                      </div>
                    ) : (
                      <div className="overflow-x-auto border border-[#E2E6EA] rounded-[6px] bg-white max-w-full">
                        <table className="w-full text-left text-[13px] border-collapse">
                          <thead>
                            <tr className="border-b border-[#E2E6EA] text-[#4B5563] bg-[#F1F3F5] text-[11px] uppercase tracking-[0.04em] select-none">
                              {[
                                { id: 'sample_id', label: 'Sample ID', align: 'left' },
                                { id: 'trimmed_reads', label: 'Filtered Reads', align: 'right' },
                                { id: 'mapped_reads', label: 'Mapped Reads', align: 'right' },
                                { id: 'mean_depth', label: 'Mean Depth', align: 'right' },
                                { id: 'coverage_breadth', label: `Cov Breadth (≥${minCovThreshold}x)`, align: 'right' },
                                { id: 'n_content', label: 'Internal N%', align: 'right' },
                                { id: 'consensus_length', label: 'Consensus Len', align: 'right' },
                                { id: 'snps_called', label: 'Variants', align: 'right' },
                              ].map((col) => {
                                const isActiveSort = sortField === col.id;
                                return (
                                  <th
                                    key={col.id}
                                    className={`${
                                      col.id === 'sample_id' ? 'pl-3.5 pr-2.5' : 'px-2.5'
                                    } py-2.5 font-semibold cursor-pointer hover:text-[#111827] transition-colors whitespace-nowrap ${
                                      col.align === 'right' ? 'text-right' : 'text-left'
                                    }`}
                                    onClick={() => handleSort(col.id)}
                                  >
                                    <div className={`flex items-center gap-1 ${col.align === 'right' ? 'justify-end' : ''}`}>
                                      {col.label}
                                      <ArrowUpDown
                                        className={`w-3 h-3 shrink-0 ${
                                          isActiveSort ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'
                                        }`}
                                      />
                                    </div>
                                  </th>
                                );
                              })}
                              <th
                                className="pl-2.5 pr-3.5 py-2.5 font-semibold cursor-pointer hover:text-[#111827] transition-colors text-center whitespace-nowrap"
                                onClick={() => handleSort('status')}
                              >
                                <div className="flex items-center justify-center gap-1">
                                  Status
                                  <ArrowUpDown
                                    className={`w-3 h-3 shrink-0 ${
                                      sortField === 'status' ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'
                                    }`}
                                  />
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
                                  className="border-b border-[#E2E6EA] last:border-none hover:bg-[#F8F9FA] cursor-pointer transition-colors"
                                  onClick={() => {
                                    setSelectedSampleId(sample.sample_id);
                                    if (refName) setSelectedRefName(refName);
                                    setActiveTab('details');
                                  }}
                                >
                                  <td className="pl-3.5 pr-2.5 py-2.5 font-semibold text-[#111827] max-w-[230px]">
                                    <div className="truncate" title={sample.sample_id}>{sample.sample_id}</div>
                                    {coDetections.length > 1 ? (
                                      <div className="text-[11px] text-[#0F5C6E] font-semibold mt-0.5 font-mono break-words leading-snug">
                                        Co-detection: {coDetections.join(' + ')}
                                      </div>
                                    ) : (
                                      refName && (
                                        <div className="text-[11px] text-[#6B7280] font-normal mt-0.5 font-mono truncate" title={refName}>
                                          {refName}
                                        </div>
                                      )
                                    )}
                                  </td>
                                  <td className="px-2.5 py-2.5 font-mono text-right text-[#111827] whitespace-nowrap">
                                    <div>{filteredReads.toLocaleString()}</div>
                                    <div className="text-[11px] text-[#6B7280]">
                                      Dup: {formatMetric(sample.qc?.duplication_rate, { decimals: 1, suffix: '%', isRatioPercent: true })}
                                    </div>
                                  </td>
                                  <td className="px-2.5 py-2.5 font-mono text-right text-[#111827] whitespace-nowrap">
                                    {mappedReads === null
                                      ? 'N/A'
                                      : `${mappedReads.toLocaleString()} (${mappedPct || '0.0%'})`}
                                  </td>
                                  <td
                                    className={`px-2.5 py-2.5 font-mono text-right font-semibold whitespace-nowrap ${
                                      evalDepth >= qcThresholds.pass_depth
                                        ? 'text-[#16654A]'
                                        : evalDepth >= qcThresholds.warn_depth
                                        ? 'text-[#8C5308]'
                                        : 'text-[#A62C22]'
                                    }`}
                                  >
                                    <div>{formatDepth(evalDepth)}</div>
                                    {cov.callable_mean_depth !== undefined && cov.mean_depth !== cov.callable_mean_depth && (
                                      <div className="text-[11px] text-[#6B7280] font-normal">
                                        whole: {formatDepth(cov.mean_depth)}
                                      </div>
                                    )}
                                  </td>
                                  <td
                                    className={`px-2.5 py-2.5 font-mono text-right font-semibold whitespace-nowrap ${
                                      evalBreadth >= qcThresholds.pass_breadth
                                        ? 'text-[#16654A]'
                                        : evalBreadth >= qcThresholds.warn_breadth
                                        ? 'text-[#8C5308]'
                                        : 'text-[#A62C22]'
                                    }`}
                                  >
                                    <div>{evalBreadth.toFixed(2)}%</div>
                                    {cov.callable_breadth_min_cov !== undefined && cov.qc_basis !== 'whole_contig' && (
                                      <div className="text-[11px] text-[#6B7280] font-normal">
                                        whole: {wholeBreadth.toFixed(2)}%
                                      </div>
                                    )}
                                  </td>
                                  <td
                                    className={`px-2.5 py-2.5 font-mono text-right font-semibold whitespace-nowrap ${
                                      evalNRate <= qcThresholds.pass_max_n
                                        ? 'text-[#16654A]'
                                        : evalNRate <= qcThresholds.warn_max_n
                                        ? 'text-[#8C5308]'
                                        : 'text-[#A62C22]'
                                    }`}
                                  >
                                    <div>{evalNRate.toFixed(2)}%</div>
                                    {cons.internal_n_content !== undefined && cov.qc_basis !== 'whole_contig' && (
                                      <div className="text-[11px] text-[#6B7280] font-normal">
                                        whole: {(cons.n_content ?? 0).toFixed(2)}%
                                      </div>
                                    )}
                                  </td>
                                  <td className="px-2.5 py-2.5 font-mono text-right text-[#111827] whitespace-nowrap">
                                    {(cons.length ?? 0).toLocaleString()} bp
                                  </td>
                                  <td className="px-2.5 py-2.5 font-mono text-right text-[#111827] whitespace-nowrap">
                                    {primaryVarsCount}
                                  </td>
                                  <td className="pl-2.5 pr-3.5 py-2.5 text-center whitespace-nowrap">
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
                    <div className="flex flex-col gap-5">
                      <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-4 p-4 glass-panel">
                        <div className="flex flex-wrap items-center gap-3">
                          <div className="flex items-center gap-2">
                            <span className="text-[11px] text-[#6B7280] font-semibold uppercase tracking-[0.06em]">
                              Active Sample:
                            </span>
                            <select
                              value={selectedSampleId}
                              onChange={(e) => {
                                const newId = e.target.value;
                                setSelectedSampleId(newId);
                                const nextSample = reportData.samples.find((s) => s.sample_id === newId);
                                const bestRef = getBestReference(nextSample, reportData);
                                if (bestRef) setSelectedRefName(bestRef.name);
                              }}
                              className="bg-white border border-[#CBD2D9] text-[#111827] text-xs font-medium rounded-[4px] px-3 py-1.5 focus:outline-none focus:border-[#0F5C6E] cursor-pointer"
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
                              <span className="text-[11px] text-[#6B7280] font-semibold uppercase tracking-[0.06em]">
                                Reference:
                              </span>
                              <select
                                value={targetRef}
                                onChange={(e) => setSelectedRefName(e.target.value)}
                                className="bg-white border border-[#CBD2D9] text-[#111827] text-xs font-mono rounded-[4px] px-3 py-1.5 focus:outline-none focus:border-[#0F5C6E] cursor-pointer"
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
                            className="flex items-center gap-1.5 px-3 py-1.5 bg-white border border-[#CBD2D9] text-[#111827] text-xs font-semibold rounded-[4px] hover:border-[#0F5C6E] hover:text-[#0F5C6E] disabled:opacity-40 transition-colors"
                          >
                            <Download className="w-3.5 h-3.5 text-[#0F5C6E]" />
                            Sample FASTA
                          </button>
                        </div>

                        <div className="flex flex-wrap items-center gap-x-3 gap-y-1 text-xs text-[#4B5563] tabular-nums">
                          <span>
                            Mapped:{' '}
                            <strong className="text-[#111827] font-mono">
                              {mappedReads === null
                                ? 'N/A'
                                : `${mappedReads.toLocaleString()} (${mappedPct || '0.0'}%)`}
                            </strong>
                          </span>
                          <span className="text-[#CBD2D9]">|</span>
                          <span>
                            Callable Depth:{' '}
                            <strong
                              className={`font-mono ${
                                evalDepth >= qcThresholds.pass_depth
                                  ? 'text-[#16654A]'
                                  : evalDepth >= qcThresholds.warn_depth
                                  ? 'text-[#8C5308]'
                                  : 'text-[#A62C22]'
                              }`}
                            >
                              {formatDepth(evalDepth)}
                            </strong>
                          </span>
                          <span className="text-[#CBD2D9]">|</span>
                          <span>
                            Breadth (&ge;{minCovThreshold}x):{' '}
                            <strong className="font-mono text-[#111827]">{evalBreadth.toFixed(2)}%</strong>
                          </span>
                          <span className="text-[#CBD2D9]">|</span>
                          <span>
                            Q30 Trim:{' '}
                            <strong className="font-mono text-[#111827]">
                              {formatMetric(activeSample.qc?.q30_rate_trimmed, { decimals: 1, suffix: '%', isRatioPercent: true })}
                            </strong>
                          </span>
                          <span className="text-[#CBD2D9]">|</span>
                          <span>
                            Dup:{' '}
                            <strong className="font-mono text-[#111827]">
                              {formatMetric(activeSample.qc?.duplication_rate, { decimals: 1, suffix: '%', isRatioPercent: true })}
                            </strong>
                          </span>
                        </div>
                      </div>

                      {activeSample.metadata && Object.keys(activeSample.metadata).length > 0 && (
                        <div className="flex flex-wrap items-center gap-x-6 gap-y-1.5 text-xs bg-[#F1F3F5] border border-[#E2E6EA] rounded-[6px] px-4 py-2.5">
                          {Object.entries(activeSample.metadata)
                            .filter(([k]) => k !== 'sample_id')
                            .map(([key, val]) => (
                              <span key={key}>
                                <span className="text-[11px] text-[#6B7280] uppercase font-semibold tracking-[0.04em] mr-1.5">
                                  {key.replace(/_/g, ' ')}:
                                </span>
                                <strong className="text-[#111827] font-semibold">{val}</strong>
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
                  const overallStatus = getSampleStatus(activeSample, reportData);
                  const stableReportId = computeStableReportId(reportData.run_id, activeSample.sample_id);

                  return (
                    <div className="flex flex-col gap-5">
                      <div className="flex flex-wrap items-center justify-between gap-4 p-4 glass-panel no-print">
                        <div>
                          <h4 className="text-[14px] font-semibold text-[#111827]">Genomic Epidemiological PDF Generator</h4>
                          <p className="text-[12px] text-[#6B7280]">
                            Multi-page print layout with QC verdict, run configuration, and software provenance.
                          </p>
                        </div>

                        <div className="flex items-center gap-3">
                          <select
                            value={selectedSampleId}
                            onChange={(e) => setSelectedSampleId(e.target.value)}
                            className="bg-white border border-[#CBD2D9] text-[#111827] text-xs font-medium rounded-[4px] px-3 py-1.5 focus:outline-none focus:border-[#0F5C6E] cursor-pointer"
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
                            className="flex items-center gap-1.5 px-3 py-1.5 bg-[#0F5C6E] hover:bg-[#0B4654] text-white text-xs font-semibold rounded-[4px] transition-colors"
                          >
                            <Printer className="w-3.5 h-3.5" />
                            Print / Save PDF
                          </button>
                        </div>
                      </div>

                      <div className="p-8 bg-white text-[#111827] rounded-[6px] max-w-4xl mx-auto flex flex-col gap-7 border border-[#E2E6EA] print-card w-full">
                        <div className="flex items-start justify-between border-b border-[#E2E6EA] pb-5 print-avoid-break">
                          <div className="flex flex-col gap-1.5">
                            <div className="flex items-center gap-2.5">
                              <span className="text-[11px] font-semibold uppercase tracking-[0.06em] text-[#0F5C6E]">
                                Genomic Epidemiology Report
                              </span>
                              <span
                                className={`inline-flex items-center gap-1 text-[11px] font-semibold px-2 py-0.5 rounded-[4px] border ${getQcBadgeClasses(
                                  overallStatus.status
                                )}`}
                              >
                                QC Verdict: {overallStatus.status}
                              </span>
                            </div>
                            <h2 className="text-[22px] font-bold text-[#111827] leading-tight">
                              Pathogen Genomics Summary
                            </h2>
                            <span className="text-xs text-[#6B7280] font-mono">Report ID: {stableReportId}</span>
                          </div>
                          <div className="flex flex-col items-end text-right gap-1 text-xs">
                            <span className="font-semibold text-[#111827]">Viromics ({platformLabel})</span>
                            <span className="text-[#4B5563] font-mono">Run ID: {reportData.run_id}</span>
                            <span className="text-[#6B7280] font-mono text-[11px]">
                              Schema v{reportData.schema_version || '1.0'} · Depth&ge;{minCovThreshold}x
                            </span>
                          </div>
                        </div>

                        <div className="grid grid-cols-2 md:grid-cols-4 gap-5 bg-[#F8F9FA] p-4 rounded-[6px] border border-[#E2E6EA] print-avoid-break">
                          <div className="flex flex-col gap-1">
                            <span className="text-[11px] uppercase font-semibold tracking-[0.04em] text-[#6B7280]">
                              Sample Identifier
                            </span>
                            <strong className="text-[13px] font-semibold text-[#111827]">{activeSample.sample_id}</strong>
                          </div>
                          <div className="flex flex-col gap-1">
                            <span className="text-[11px] uppercase font-semibold tracking-[0.04em] text-[#6B7280]">
                              Reference Target
                            </span>
                            <strong className="text-[13px] font-semibold text-[#111827] font-mono break-all">
                              {coDetections.length > 1 ? `Co-detection: ${coDetections.join(' + ')}` : pdfRefName}
                            </strong>
                          </div>
                          <div className="flex flex-col gap-1">
                            <span className="text-[11px] uppercase font-semibold tracking-[0.04em] text-[#6B7280]">
                              Consensus Length
                            </span>
                            <strong className="text-[13px] font-semibold text-[#111827] font-mono tabular-nums">
                              {coDetections.length > 1
                                ? coDetections
                                    .map((c) => `${(getContigConsensus(activeSample, c).length ?? 0).toLocaleString()} bp (${c})`)
                                    .join(' | ')
                                : `${(pdfCons.length ?? 0).toLocaleString()} bp`}
                            </strong>
                          </div>
                          <div className="flex flex-col gap-1">
                            <span className="text-[11px] uppercase font-semibold tracking-[0.04em] text-[#6B7280]">
                              Evaluation Basis
                            </span>
                            <strong className="text-[13px] font-semibold text-[#111827] font-mono">
                              {overallStatus.basis || 'whole_contig'}
                            </strong>
                          </div>
                        </div>

                        <div className="print-avoid-break">
                          <h3 className="text-[11px] font-semibold text-[#111827] border-b border-[#E2E6EA] pb-2 mb-3.5 uppercase tracking-[0.06em]">
                            Sequencing Quality Assurance
                          </h3>
                          {coDetections.length > 1 ? (
                            <div className="overflow-x-auto border border-[#E2E6EA] rounded-[6px]">
                              <table className="w-full text-left text-xs border-collapse">
                                <thead>
                                  <tr className="bg-[#F1F3F5] text-[#4B5563] font-semibold border-b border-[#E2E6EA] text-[11px] uppercase tracking-[0.04em]">
                                    <th className="p-3">Detected Strain</th>
                                    <th className="p-3 text-right">Mapped Reads</th>
                                    <th className="p-3 text-right">Callable Depth</th>
                                    <th className="p-3 text-right">Callable Breadth (&ge;{minCovThreshold}x)</th>
                                    <th className="p-3 text-right">Internal N Rate</th>
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
                                      <tr key={c} className="border-b border-[#E2E6EA] last:border-none text-[#111827]">
                                        <td className="p-3 font-mono font-semibold text-[#111827]">{c}</td>
                                        <td className="p-3 font-mono text-right whitespace-nowrap">
                                          {mr === null ? 'N/A' : `${mr.toLocaleString()} (${pct || '0.0'}%)`}
                                        </td>
                                        <td className="p-3 font-mono text-right font-semibold whitespace-nowrap">
                                          {formatDepth(depth)}
                                        </td>
                                        <td className="p-3 font-mono text-right font-semibold whitespace-nowrap">
                                          {breadth.toFixed(2)}%
                                        </td>
                                        <td className="p-3 font-mono text-right font-semibold whitespace-nowrap">
                                          {nRate.toFixed(2)}%
                                        </td>
                                        <td className="p-3 text-center whitespace-nowrap">
                                          <span
                                            className={`inline-block px-2 py-0.5 rounded-[4px] border text-[11px] font-semibold ${getQcBadgeClasses(
                                              statusInfo.status
                                            )}`}
                                          >
                                            {statusInfo.status}
                                          </span>
                                        </td>
                                      </tr>
                                    );
                                  })}
                                </tbody>
                              </table>
                            </div>
                          ) : (
                            <div className="grid grid-cols-2 md:grid-cols-5 gap-3 text-center">
                              <div className="p-3 border border-[#E2E6EA] rounded-[4px] bg-[#F8F9FA]">
                                <span className="text-[10px] uppercase font-semibold tracking-[0.04em] text-[#6B7280] block">
                                  Filtered Reads
                                </span>
                                <strong className="text-sm text-[#111827] font-semibold font-mono tabular-nums mt-0.5 block">
                                  {(activeSample.qc?.trimmed_reads ?? 0).toLocaleString()}
                                </strong>
                              </div>
                              <div className="p-3 border border-[#E2E6EA] rounded-[4px] bg-[#F8F9FA]">
                                <span className="text-[10px] uppercase font-semibold tracking-[0.04em] text-[#6B7280] block">
                                  Mapped Reads
                                </span>
                                <strong className="text-sm text-[#111827] font-semibold font-mono tabular-nums mt-0.5 block">
                                  {(() => {
                                    const mr = getMappedReads(activeSample, pdfRefName);
                                    return mr === null ? 'N/A' : mr.toLocaleString();
                                  })()}
                                </strong>
                              </div>
                              <div className="p-3 border border-[#E2E6EA] rounded-[4px] bg-[#F8F9FA]">
                                <span className="text-[10px] uppercase font-semibold tracking-[0.04em] text-[#6B7280] block">
                                  Callable Depth
                                </span>
                                <strong className="text-sm text-[#111827] font-semibold font-mono tabular-nums mt-0.5 block">
                                  {formatDepth(getEvalDepth(activeSample, pdfRefName))}
                                </strong>
                              </div>
                              <div className="p-3 border border-[#E2E6EA] rounded-[4px] bg-[#F8F9FA]">
                                <span className="text-[10px] uppercase font-semibold tracking-[0.04em] text-[#6B7280] block">
                                  Breadth (&ge;{minCovThreshold}x)
                                </span>
                                <strong className="text-sm text-[#111827] font-semibold font-mono tabular-nums mt-0.5 block">
                                  {getEvalBreadth(activeSample, pdfRefName).toFixed(2)}%
                                </strong>
                              </div>
                              <div className="p-3 border border-[#E2E6EA] rounded-[4px] bg-[#F8F9FA]">
                                <span className="text-[10px] uppercase font-semibold tracking-[0.04em] text-[#6B7280] block">
                                  Internal N Rate
                                </span>
                                <strong className="text-sm text-[#111827] font-semibold font-mono tabular-nums mt-0.5 block">
                                  {getEvalNRate(activeSample, pdfRefName).toFixed(2)}%
                                </strong>
                                {pdfCov.qc_basis && pdfCov.qc_basis !== 'whole_contig' && (
                                  <span className="text-[10px] text-[#6B7280] font-mono block">
                                    whole N: {(pdfCons.n_content ?? 0).toFixed(2)}%
                                  </span>
                                )}
                              </div>
                            </div>
                          )}
                        </div>

                        {versionEntries.length > 0 && (
                          <div className="bg-[#F1F3F5] border border-[#E2E6EA] rounded-[4px] p-3 text-[11px] font-mono text-[#4B5563] print-avoid-break">
                            <strong className="text-[#111827] mr-2">Provenance:</strong>
                            {versionEntries.map(([tool, ver]) => `${tool}=${ver}`).join(' · ')}
                          </div>
                        )}

                        <div className="flex flex-col gap-5">
                          {(() => {
                            const allVariants = activeSample.variants || [];
                            const renderMutationsTable = (titleLabel, filteredVars) => (
                              <div className="flex flex-col gap-3">
                                <div className="flex justify-between items-center border-b border-[#E2E6EA] pb-2">
                                  <h3 className="text-[11px] font-semibold text-[#111827]">
                                    <span className="uppercase tracking-[0.06em]">Identified Variants:</span>{' '}
                                    <span className="font-mono normal-case font-semibold text-[#0F5C6E]">
                                      {titleLabel}
                                    </span>
                                  </h3>
                                  <span className="text-[11px] text-[#4B5563] font-mono font-semibold bg-[#F1F3F5] border border-[#E2E6EA] px-2 py-0.5 rounded-[4px] tabular-nums">
                                    {filteredVars.length} variants
                                  </span>
                                </div>
                                <div className="overflow-x-auto border border-[#E2E6EA] rounded-[6px]">
                                  <table className="w-full text-left text-xs border-collapse">
                                    <thead>
                                      <tr className="bg-[#F1F3F5] border-b border-[#E2E6EA] text-[#4B5563] font-semibold text-[11px] uppercase tracking-[0.04em]">
                                        <th className="p-3 whitespace-nowrap">Position</th>
                                        <th className="p-3 whitespace-nowrap">Allele</th>
                                        <th className="p-3 whitespace-nowrap">Gene Region</th>
                                        <th className="p-3">Protein Change</th>
                                        <th className="p-3 whitespace-nowrap">Type</th>
                                        <th className="p-3 text-right whitespace-nowrap">Allele Freq</th>
                                      </tr>
                                    </thead>
                                    <tbody>
                                      {filteredVars.length > 0 ? (
                                        filteredVars.map((v, i) => {
                                          const cons = deduplicateConsequences(formatConsequences(v));
                                          return (
                                            <tr key={i} className="border-b border-[#E2E6EA] last:border-none text-[#111827]">
                                              <td className="p-3 font-mono font-medium whitespace-nowrap">
                                                {v.position.toLocaleString()}
                                              </td>
                                              <td className="p-3 font-mono font-semibold text-[#111827] whitespace-nowrap">
                                                {v.ref} &gt; {v.alt}
                                              </td>
                                              <td className="p-3 whitespace-nowrap">
                                                <span className="px-2 py-0.5 rounded-[4px] text-[11px] font-semibold bg-[#F1F3F5] text-[#111827] border border-[#E2E6EA]">
                                                  {v.gene}
                                                </span>
                                              </td>
                                              <td className="p-3 font-mono text-[11px]">
                                                {cons.length > 0 ? (
                                                  <div className="flex flex-wrap gap-1.5">
                                                    {cons.map((c, cIdx) => (
                                                      <span
                                                        key={cIdx}
                                                        title={c.tooltip}
                                                        className={
                                                          c.kind === 'unavailable'
                                                            ? 'text-[#6B7280] italic'
                                                            : c.kind === 'synonymous'
                                                            ? 'text-[#5C6773]'
                                                            : c.kind === 'frameshift' || c.kind === 'indel'
                                                            ? 'text-[#9A5B06] font-semibold'
                                                            : 'text-[#B43B2E] font-semibold'
                                                        }
                                                      >
                                                        {cons.length > 1 ? `${c.gene}: ` : ''}
                                                        {c.label}
                                                      </span>
                                                    ))}
                                                  </div>
                                                ) : (
                                                  <span className="text-[#6B7280]">-</span>
                                                )}
                                              </td>
                                              <td className="p-3 whitespace-nowrap text-[#4B5563]">{v.type}</td>
                                              <td className="p-3 font-mono text-right text-[#0F5C6E] font-semibold whitespace-nowrap">
                                                {(v.frequency * 100).toFixed(1)}%
                                              </td>
                                            </tr>
                                          );
                                        })
                                      ) : (
                                        <tr>
                                          <td colSpan={6} className="p-4 text-center text-[#6B7280]">
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

                        <div className="flex justify-between items-center border-t border-[#E2E6EA] pt-7 mt-2 text-xs print-avoid-break">
                          <div className="flex flex-col gap-1">
                            <div className="w-36 border-b border-[#CBD2D9] h-8" />
                            <span className="text-[#6B7280] uppercase font-semibold text-[10px] tracking-[0.06em]">
                              Clinician / Assayer Signature
                            </span>
                          </div>
                          <div className="flex flex-col gap-1 text-right">
                            <div className="w-36 border-b border-[#CBD2D9] h-8 self-end" />
                            <span className="text-[#6B7280] uppercase font-semibold text-[10px] tracking-[0.06em]">
                              Lab Director Approval
                            </span>
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
