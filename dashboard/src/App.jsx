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
  FileCode,
  ShieldCheck,
  Zap,
  Printer,
  ArrowUpDown
} from 'lucide-react';
import confetti from 'canvas-confetti';

// Import custom components
import CoverageChart from './components/CoverageChart';
import GeneMap from './components/GeneMap';
import MutationMatrix from './components/MutationMatrix';
import VariantTable from './components/VariantTable';
import PhylogeneticTree from './components/PhylogeneticTree';

export default function App() {
  const [reportData, setReportData] = useState(null);
  const [selectedSampleId, setSelectedSampleId] = useState('');
  const [activeTab, setActiveTab] = useState('overview'); // 'overview', 'details', 'matrix', 'tree', 'report'
  const [isDragging, setIsDragging] = useState(false);
  const [uploadError, setUploadError] = useState('');
  const [sortField, setSortField] = useState('sample_id');
  const [sortAsc, setSortAsc] = useState(true);
  const [selectedRefName, setSelectedRefName] = useState('');

  const minCovThreshold = reportData ? (reportData.min_coverage_threshold || 10) : 10;

  // Load mock data on mount as the interactive demo
  useEffect(() => {
    fetch('/mock_run_report.json')
      .then(res => {
        if (!res.ok) throw new Error('Demo report not found');
        return res.json();
      })
      .then(data => {
        setReportData(data);
        if (data.samples && data.samples.length > 0) {
          setSelectedSampleId(data.samples[0].sample_id);
        }
        if (data.references && data.references.length > 0) {
          setSelectedRefName(data.references[0].name);
        } else if (data.reference) {
          setSelectedRefName(data.reference.name);
        }
      })
      .catch(err => {
        console.warn('Could not load demo report automatically:', err);
      });
  }, []);

  // Trigger celebration on successful custom file upload
  const celebrateSuccess = () => {
    confetti({
      particleCount: 80,
      spread: 60,
      origin: { y: 0.8 },
      colors: ['#14b8a6', '#8b5cf6', '#10b981']
    });
  };

  // File Upload Handlers
  const validateUploadJSON = (data) => {
    if (!data || typeof data !== 'object') {
      throw new Error('Report must be a JSON object.');
    }
    if (!data.run_id || typeof data.run_id !== 'string') {
      throw new Error("Missing or invalid 'run_id' in report.");
    }
    if (data.min_coverage_threshold === undefined || typeof data.min_coverage_threshold !== 'number') {
      throw new Error("Missing or invalid 'min_coverage_threshold' in report.");
    }
    if (!data.references || !Array.isArray(data.references)) {
      throw new Error("Missing or invalid 'references' list in report.");
    }
    if (!data.samples || !Array.isArray(data.samples)) {
      throw new Error("Missing or invalid 'samples' list in report.");
    }
    
    // Check references structure
    data.references.forEach((ref, idx) => {
      if (!ref.name || typeof ref.name !== 'string') {
        throw new Error(`Reference at index ${idx} is missing 'name'.`);
      }
      if (!ref.length || typeof ref.length !== 'number') {
        throw new Error(`Reference '${ref.name}' is missing 'length'.`);
      }
    });

    // Check samples structure
    data.samples.forEach((s, idx) => {
      const sId = s.sample_id || `index_${idx}`;
      if (!s.sample_id || typeof s.sample_id !== 'string') {
        throw new Error(`Sample at index ${idx} is missing 'sample_id'.`);
      }
      if (!s.qc || typeof s.qc !== 'object') {
        throw new Error(`Sample '${sId}' is missing 'qc' dictionary.`);
      }
      if (!s.coverage || typeof s.coverage !== 'object') {
        throw new Error(`Sample '${sId}' is missing 'coverage' dictionary.`);
      }
      if (!s.consensus || typeof s.consensus !== 'object') {
        throw new Error(`Sample '${sId}' is missing 'consensus' dictionary.`);
      }
      
      // Check qc sub-keys
      const requiredQc = ["raw_reads", "raw_bases", "trimmed_reads", "trimmed_bases", "q30_rate_raw", "q30_rate_trimmed", "duplication_rate"];
      requiredQc.forEach(key => {
        if (s.qc[key] === undefined || typeof s.qc[key] !== 'number') {
          throw new Error(`Sample '${sId}' QC is missing or has invalid number for '${key}'.`);
        }
      });

      // Check coverage sub-keys
      Object.keys(s.coverage).forEach(refName => {
        const cov = s.coverage[refName];
        if (!cov || typeof cov !== 'object') {
          throw new Error(`Sample '${sId}' is missing coverage data for reference '${refName}'.`);
        }
        const requiredCov = ["mean_depth", "coverage_breadth_10x", "coverage_breadth_30x", "mapped_reads", "coverage_profile"];
        requiredCov.forEach(key => {
          if (cov[key] === undefined) {
            throw new Error(`Sample '${sId}' coverage for reference '${refName}' is missing '${key}'.`);
          }
        });
        if (!Array.isArray(cov.coverage_profile)) {
          throw new Error(`Sample '${sId}' coverage profile for reference '${refName}' must be a list.`);
        }
      });

      // Check consensus sub-keys
      Object.keys(s.consensus).forEach(refName => {
        const cons = s.consensus[refName];
        if (!cons || typeof cons !== 'object') {
          throw new Error(`Sample '${sId}' is missing consensus data for reference '${refName}'.`);
        }
        const requiredCons = ["length", "gc_content", "n_content", "sequence"];
        requiredCons.forEach(key => {
          if (cons[key] === undefined) {
            throw new Error(`Sample '${sId}' consensus for reference '${refName}' is missing '${key}'.`);
          }
        });
      });
    });
  };

  const handleJSONData = (text) => {
    try {
      const data = JSON.parse(text);
      validateUploadJSON(data);
      setReportData(data);
      if (data.samples.length > 0) {
        setSelectedSampleId(data.samples[0].sample_id);
      }
      if (data.references && data.references.length > 0) {
        setSelectedRefName(data.references[0].name);
      } else if (data.reference) {
        setSelectedRefName(data.reference.name);
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
    if (file && file.type === "application/json") {
      const reader = new FileReader();
      reader.onload = (event) => handleJSONData(event.target.result);
      reader.readAsText(file);
    } else {
      setUploadError('Only JSON run report files are supported.');
    }
  };

  if (!reportData) {
    return (
      <div className="flex flex-col items-center justify-center min-h-screen p-6">
        <div className="w-full max-w-lg p-8 glass-panel rounded-3xl flex flex-col gap-6 text-center shadow-glow-teal animate-pulse">
          <Dna className="w-16 h-16 mx-auto text-brand-teal animate-bounce" />
          <h2 className="text-2xl font-bold font-outfit text-white">Loading Viromics Suite...</h2>
          <p className="text-xs text-gray-400">Initializing analytical frameworks and visuals.</p>
        </div>
      </div>
    );
  }

  // Active sample computed properties
  const activeSample = reportData.samples.find(s => s.sample_id === selectedSampleId) || reportData.samples[0];

  // Get all references present in the report
  const getReferencesList = () => {
    if (reportData.references) return reportData.references;
    if (reportData.reference) return [reportData.reference];
    return [];
  };

  // Helper: check sample status for a specific reference
  const getSampleStatusForRef = (sample, refName) => {
    if (!sample.coverage) return 'FAIL';
    
    let cov = sample.coverage;
    if (cov[refName]) {
      cov = cov[refName];
    }
    
    const breadth = cov.coverage_breadth_10x !== undefined ? cov.coverage_breadth_10x : 0;
    const reads = sample.qc.trimmed_reads;
    
    if (breadth >= 95 && reads >= 5000) return 'PASS';
    if (breadth >= 80 && reads >= 2000) return 'WARN';
    return 'FAIL';
  };

  // Helper: Get the best matching reference for a sample (highest coverage breadth)
  const getBestReference = (sample) => {
    const refs = getReferencesList();
    if (refs.length === 0) return null;
    if (refs.length === 1) return refs[0];
    
    let bestRef = refs[0];
    let maxBreadth = -1;
    
    refs.forEach(r => {
      const cov = sample.coverage && sample.coverage[r.name];
      const breadth = cov ? cov.coverage_breadth_10x : 0;
      if (breadth > maxBreadth) {
        maxBreadth = breadth;
        bestRef = r;
      }
    });
    
    return bestRef;
  };

  // Helper: Get the overall status of the sample (based on its best reference)
  const getSampleStatus = (sample) => {
    const bestRef = getBestReference(sample);
    if (!bestRef) return 'FAIL';
    return getSampleStatusForRef(sample, bestRef.name);
  };

  // Helper: Detect co-infections/co-detections
  // Defined as having more than one strain with >= 10% coverage breadth
  const getCoDetections = (sample) => {
    const refs = getReferencesList();
    if (refs.length <= 1) return [];
    
    const detected = [];
    refs.forEach(r => {
      const cov = sample.coverage && sample.coverage[r.name];
      const breadth = cov ? cov.coverage_breadth_10x : 0;
      if (breadth >= 10.0) {
        detected.push(r.name);
      }
    });
    return detected;
  };
  
  // Helper: Get the number of mapped reads for a sample and specific reference strain
  const getMappedReadsForSampleRef = (sample, refName) => {
    if (!sample) return 0;
    const cov = (refName && sample.coverage && sample.coverage[refName]) ? sample.coverage[refName] : (sample.coverage || {});
    if (cov.mapped_reads !== undefined) return cov.mapped_reads;
    
    if (!cov.mean_depth) return 0;
    
    const refs = getReferencesList();
    const ref = (refName && refs.find(r => r.name === refName)) || refs[0] || reportData.reference;
    const refLength = ref ? ref.length : 0;
    if (!refLength) return 0;
    
    // Estimate mapped reads: (mean_depth * refLength) / average_read_length (assumed 150bp)
    const trimmed = sample.qc?.trimmed_reads || 100000;
    const est = Math.round((cov.mean_depth * refLength) / 150);
    // Cap at a reasonable percentage of trimmed reads (e.g. 98.5%) or at least match the proportion of coverage
    const breadth = cov.coverage_breadth_10x !== undefined ? cov.coverage_breadth_10x / 100 : 0.95;
    const cappedMax = Math.round(trimmed * Math.min(0.995, breadth + 0.05));
    return Math.max(1, Math.min(est, cappedMax));
  };

  const handleSort = (field) => {
    if (sortField === field) {
      setSortAsc(!sortAsc);
    } else {
      setSortField(field);
      setSortAsc(true);
    }
  };

  // Sort samples list based on sortField (with multi-reference support)
  const sortedSamples = [...reportData.samples].sort((a, b) => {
    if (!sortField) return 0;
    
    let valA, valB;
    if (sortField === 'sample_id') {
      valA = a.sample_id;
      valB = b.sample_id;
    } else if (sortField === 'trimmed_reads') {
      valA = a.qc.trimmed_reads;
      valB = b.qc.trimmed_reads;
    } else if (sortField === 'mapped_reads') {
      const refA = getBestReference(a);
      const refB = getBestReference(b);
      valA = getMappedReadsForSampleRef(a, refA?.name);
      valB = getMappedReadsForSampleRef(b, refB?.name);
    } else if (sortField === 'mean_depth') {
      const refA = getBestReference(a);
      const refB = getBestReference(b);
      const covA = refA && a.coverage ? (a.coverage[refA.name] ? a.coverage[refA.name].mean_depth : a.coverage.mean_depth) : 0;
      const covB = refB && b.coverage ? (b.coverage[refB.name] ? b.coverage[refB.name].mean_depth : b.coverage.mean_depth) : 0;
      valA = covA || 0;
      valB = covB || 0;
    } else if (sortField === 'coverage_breadth_10x') {
      const refA = getBestReference(a);
      const refB = getBestReference(b);
      const covA = refA && a.coverage ? (a.coverage[refA.name] ? a.coverage[refA.name].coverage_breadth_10x : a.coverage.coverage_breadth_10x) : 0;
      const covB = refB && b.coverage ? (b.coverage[refB.name] ? b.coverage[refB.name].coverage_breadth_10x : b.coverage.coverage_breadth_10x) : 0;
      valA = covA || 0;
      valB = covB || 0;
    } else if (sortField === 'consensus_length') {
      const refA = getBestReference(a);
      const refB = getBestReference(b);
      const consA = refA && a.consensus ? (a.consensus[refA.name] ? a.consensus[refA.name].length : a.consensus.length) : 0;
      const consB = refB && b.consensus ? (b.consensus[refB.name] ? b.consensus[refB.name].length : b.consensus.length) : 0;
      valA = consA || 0;
      valB = consB || 0;
    } else if (sortField === 'n_content') {
      const refA = getBestReference(a);
      const refB = getBestReference(b);
      const consA = refA && a.consensus ? (a.consensus[refA.name] ? a.consensus[refA.name].n_content : a.consensus.n_content) : 0;
      const consB = refB && b.consensus ? (b.consensus[refB.name] ? b.consensus[refB.name].n_content : b.consensus.n_content) : 0;
      valA = consA || 0;
      valB = consB || 0;
    } else if (sortField === 'snps_called') {
      valA = a.variants ? a.variants.length : 0;
      valB = b.variants ? b.variants.length : 0;
    } else if (sortField === 'status') {
      valA = getSampleStatus(a);
      valB = getSampleStatus(b);
    }
    
    if (typeof valA === 'string') {
      valA = valA.toLowerCase();
      valB = valB.toLowerCase();
    }
    
    if (valA < valB) return sortAsc ? -1 : 1;
    if (valA > valB) return sortAsc ? 1 : -1;
    return 0;
  });

  return (
    <div className="min-h-screen pb-16 flex flex-col">
      
      {/* Top sleeks headers */}
      <header className="no-print border-b border-dark-border/60 bg-dark-bg/40 backdrop-blur-md sticky top-0 z-50 transition-all duration-300">
        <div className="max-w-7xl mx-auto px-6 h-16 flex items-center justify-between">
          <div className="flex items-center gap-3">
            <div className="w-9 h-9 rounded-xl bg-gradient-to-tr from-brand-teal to-brand-violet flex items-center justify-center shadow-glow-teal">
              <Dna className="w-5 h-5 text-dark-bg font-bold" />
            </div>
            <div>
              <h1 className="text-md font-bold font-outfit text-white tracking-wide flex items-center gap-1.5">
                VIROMICS <span className="text-[10px] py-0.5 px-1.5 rounded bg-brand-teal/20 text-brand-teal font-semibold font-mono tracking-normal">v1.0</span>
              </h1>
              <p className="text-[10px] text-gray-400">Genomic Epidemiology Suite</p>
            </div>
          </div>
          
          <div className="flex items-center gap-4">
            {/* Run summary pill */}
            <div className="hidden md:flex items-center gap-2 bg-dark-card border border-dark-border px-3 py-1.5 rounded-xl text-xs">
              <Zap className="w-3.5 h-3.5 text-brand-amber animate-pulse" />
              <span className="text-gray-400">Active Run:</span>
              <strong className="text-white font-mono">{reportData.run_id}</strong>
            </div>

            {/* Custom file upload button */}
            <label className="flex items-center gap-1.5 px-3 py-1.5 bg-dark-card border border-dark-border text-xs text-gray-300 font-semibold rounded-xl hover:bg-dark-hover cursor-pointer active:scale-95 transition-all">
              <UploadCloud className="w-4 h-4 text-brand-teal" />
              Upload Run JSON
              <input type="file" accept=".json" onChange={handleFileUpload} className="hidden" />
            </label>
          </div>
        </div>
      </header>

      {/* Main Container */}
      <main className="max-w-7xl mx-auto px-6 w-full flex-1 flex flex-col gap-6 mt-8">
        
        {/* Drag and Drop Zone Alert if dragging */}
        {isDragging ? (
          <div 
            onDragOver={handleDragOver}
            onDragLeave={handleDragLeave}
            onDrop={handleDrop}
            className="border-2 border-dashed border-brand-teal bg-brand-teal/10 rounded-3xl p-16 flex flex-col items-center justify-center text-center gap-4 transition-all"
          >
            <UploadCloud className="w-16 h-16 text-brand-teal animate-bounce" />
            <h2 className="text-xl font-bold text-white">Drop the report JSON file here</h2>
            <p className="text-sm text-gray-400">Instantly populate coverage, mutations and clustering profiles.</p>
          </div>
        ) : (
          <div 
            onDragOver={handleDragOver}
            className="flex flex-col gap-6"
          >
            {/* Display custom upload error */}
            {uploadError && (
              <div className="bg-brand-rose/10 border border-brand-rose/20 text-brand-rose text-xs p-3 rounded-xl flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 shrink-0" />
                <span>{uploadError}</span>
              </div>
            )}

            {/* Top Navigation Tabs */}
            <div className="flex flex-wrap border-b border-dark-border/40 gap-1 sm:gap-2 no-print">
              <button
                onClick={() => setActiveTab('overview')}
                className={`flex items-center gap-2 px-4 py-2.5 text-xs font-semibold rounded-t-xl transition-all border-b-2 ${
                  activeTab === 'overview'
                    ? 'border-brand-teal text-brand-teal bg-brand-teal/5'
                    : 'border-transparent text-gray-400 hover:text-white'
                }`}
              >
                <Grid className="w-4 h-4" />
                Run Overview
              </button>
              <button
                onClick={() => setActiveTab('details')}
                className={`flex items-center gap-2 px-4 py-2.5 text-xs font-semibold rounded-t-xl transition-all border-b-2 ${
                  activeTab === 'details'
                    ? 'border-brand-teal text-brand-teal bg-brand-teal/5'
                    : 'border-transparent text-gray-400 hover:text-white'
                }`}
              >
                <Sliders className="w-4 h-4" />
                Sample Details
              </button>
              <button
                onClick={() => setActiveTab('matrix')}
                className={`flex items-center gap-2 px-4 py-2.5 text-xs font-semibold rounded-t-xl transition-all border-b-2 ${
                  activeTab === 'matrix'
                    ? 'border-brand-teal text-brand-teal bg-brand-teal/5'
                    : 'border-transparent text-gray-400 hover:text-white'
                }`}
              >
                <Layers className="w-4 h-4" />
                Mutational Sharing
              </button>
              <button
                onClick={() => setActiveTab('tree')}
                className={`flex items-center gap-2 px-4 py-2.5 text-xs font-semibold rounded-t-xl transition-all border-b-2 ${
                  activeTab === 'tree'
                    ? 'border-brand-teal text-brand-teal bg-brand-teal/5'
                    : 'border-transparent text-gray-400 hover:text-white'
                }`}
              >
                <GitBranch className="w-4 h-4" />
                Phylogenetic Tree
              </button>
              <button
                onClick={() => setActiveTab('report')}
                className={`flex items-center gap-2 px-4 py-2.5 text-xs font-semibold rounded-t-xl transition-all border-b-2 ${
                  activeTab === 'report'
                    ? 'border-brand-teal text-brand-teal bg-brand-teal/5'
                    : 'border-transparent text-gray-400 hover:text-white'
                }`}
              >
                <FileText className="w-4 h-4" />
                Genomics Report
              </button>
            </div>

            {/* TAB CONTENT: RUN OVERVIEW */}
            {activeTab === 'overview' && (
              <div className="flex flex-col gap-6">
                
                {/* Aggregate Run cards */}
                <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4">
                  
                  {/* Card 1: Total Samples */}
                  <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                    <div>
                      <span className="text-[10px] text-gray-400 font-bold uppercase tracking-wider font-outfit">Total Samples</span>
                      <h3 className="text-2xl font-bold text-white mt-1 font-outfit">{reportData.samples.length}</h3>
                      <span className="text-[10px] text-gray-400 mt-1 block">Involved in sequencer run</span>
                    </div>
                    <div className="w-10 h-10 rounded-xl bg-brand-teal/10 flex items-center justify-center">
                      <BarChart3 className="w-5 h-5 text-brand-teal" />
                    </div>
                  </div>

                  {/* Card 2: Passing rate */}
                  {(() => {
                    const passing = reportData.samples.filter(s => getSampleStatus(s) === 'PASS').length;
                    const warn = reportData.samples.filter(s => getSampleStatus(s) === 'WARN').length;
                    const fail = reportData.samples.filter(s => getSampleStatus(s) === 'FAIL').length;
                    
                    return (
                      <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                        <div>
                          <span className="text-[10px] text-gray-400 font-bold uppercase tracking-wider font-outfit">Quality Pass Rate</span>
                          <h3 className="text-2xl font-bold text-brand-emerald mt-1 font-outfit">
                            {Math.round((passing / reportData.samples.length) * 100)}%
                          </h3>
                          <span className="text-[10px] text-gray-400 mt-1 block">
                            {passing} PASS, {warn} WARN, {fail} FAIL
                          </span>
                        </div>
                        <div className="w-10 h-10 rounded-xl bg-brand-emerald/10 flex items-center justify-center">
                          <ShieldCheck className="w-5 h-5 text-brand-emerald" />
                        </div>
                      </div>
                    );
                  })()}

                  {/* Card 3: Reference Target */}
                  <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                    <div>
                      <span className="text-[10px] text-gray-400 font-bold uppercase tracking-wider font-outfit">Reference Target</span>
                      <h3 className="text-sm font-bold text-white mt-1 font-mono truncate max-w-[160px]">
                        {reportData.reference ? reportData.reference.name : (reportData.references && reportData.references.length > 0 ? (reportData.references.length === 1 ? reportData.references[0].name : `${reportData.references.length} Reference Strains`) : 'Unknown')}
                      </h3>
                      <span className="text-[10px] text-gray-400 mt-1.5 block">
                        Length: {reportData.reference ? `${reportData.reference.length.toLocaleString()} bp` : (reportData.references && reportData.references.length > 0 ? `${reportData.references.map(r => r.length).reduce((a, b) => a + b, 0).toLocaleString()} bp total` : '0 bp')}
                      </span>
                    </div>
                    <div className="w-10 h-10 rounded-xl bg-brand-indigo/10 flex items-center justify-center">
                      <Dna className="w-5 h-5 text-brand-indigo" />
                    </div>
                  </div>

                  {/* Card 4: Run Mode */}
                  <div className="p-5 glass-panel rounded-2xl flex items-center justify-between">
                    <div>
                      <span className="text-[10px] text-gray-400 font-bold uppercase tracking-wider font-outfit font-medium">Pipeline Mode</span>
                      <h3 className="text-md font-bold text-brand-teal mt-1.5 font-outfit uppercase">
                        Amplicon WGS
                      </h3>
                      <span className="text-[10px] text-gray-400 mt-2 block">
                        Illumina + ONT support
                      </span>
                    </div>
                    <div className="w-10 h-10 rounded-xl bg-brand-teal/10 flex items-center justify-center">
                      <Activity className="w-5 h-5 text-brand-teal animate-pulse" />
                    </div>
                  </div>
                </div>

                {/* Samples comparison list */}
                <div className="p-6 glass-panel rounded-2xl flex flex-col gap-4">
                  <div>
                    <h3 className="text-lg font-semibold font-outfit text-white">Consensus Summary & QC</h3>
                    <p className="text-xs text-gray-400">
                      High-level run comparison metrics across all demultiplexed files
                    </p>
                  </div>
                  
                  <div className="overflow-x-auto border border-dark-border rounded-xl bg-dark-bg/40 max-w-full">
                    <table className="w-full text-left text-xs border-collapse">
                      <thead>
                        <tr className="border-b border-dark-border text-gray-400 bg-dark-bg/60 select-none">
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('sample_id')}>
                            <div className="flex items-center gap-1">
                              Sample ID <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('trimmed_reads')}>
                            <div className="flex items-center gap-1">
                              Filtered Reads <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('mapped_reads')}>
                            <div className="flex items-center gap-1">
                              Mapped Reads <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('mean_depth')}>
                            <div className="flex items-center gap-1">
                              Mean Depth <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('coverage_breadth_10x')}>
                            <div className="flex items-center gap-1">
                              Cov Breadth (≥{minCovThreshold}x) <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('consensus_length')}>
                            <div className="flex items-center gap-1">
                              Consensus Length <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('n_content')}>
                            <div className="flex items-center gap-1">
                              N Content <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors" onClick={() => handleSort('snps_called')}>
                            <div className="flex items-center gap-1">
                              SNPs called <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                          <th className="p-4 font-semibold font-outfit cursor-pointer hover:text-white transition-colors text-center" onClick={() => handleSort('status')}>
                            <div className="flex items-center justify-center gap-1">
                              Status <ArrowUpDown className="w-3 h-3 text-brand-teal/70" />
                            </div>
                          </th>
                        </tr>
                      </thead>
                      <tbody>
                        {sortedSamples.map((sample) => {
                          const status = getSampleStatus(sample);
                          const bestRef = getBestReference(sample);
                          const cov = (bestRef && sample.coverage && sample.coverage[bestRef.name]) ? sample.coverage[bestRef.name] : (sample.coverage || {});
                          const cons = (bestRef && sample.consensus && sample.consensus[bestRef.name]) ? sample.consensus[bestRef.name] : (sample.consensus || {});
                          const coDetections = getCoDetections(sample);
                          const totalSNPs = sample.variants ? sample.variants.length : 0;
                          
                          return (
                            <tr 
                              key={sample.sample_id} 
                              className="border-b border-dark-border/40 hover:bg-dark-hover/40 cursor-pointer transition-colors"
                              onClick={() => {
                                setSelectedSampleId(sample.sample_id);
                                setActiveTab('details');
                              }}
                            >
                              <td className="p-4 font-semibold text-white font-outfit">
                                <div>{sample.sample_id}</div>
                                {coDetections.length > 1 ? (
                                  <div className="text-[10px] text-brand-violet font-bold mt-0.5">
                                    Co-detection: {coDetections.join(' + ')}
                                  </div>
                                ) : (
                                  bestRef && (
                                    <div className="text-[10px] text-gray-400 font-normal mt-0.5">
                                      Strain: {bestRef.name}
                                    </div>
                                  )
                                )}
                              </td>
                              <td className="p-4 font-mono text-gray-300">
                                {sample.qc.trimmed_reads.toLocaleString()}
                              </td>
                              <td className="p-4 font-mono text-brand-indigo/90 font-medium">
                                {(() => {
                                  const mr = getMappedReadsForSampleRef(sample, bestRef?.name);
                                  const pct = sample.qc.trimmed_reads > 0 ? ((mr / sample.qc.trimmed_reads) * 100).toFixed(1) : '0.0';
                                  return `${mr.toLocaleString()} (${pct}%)`;
                                })()}
                              </td>
                              <td className="p-4 font-mono text-gray-300 font-semibold">
                                {cov.mean_depth !== undefined ? cov.mean_depth.toLocaleString() : 0}x
                              </td>
                              <td className="p-4 font-mono text-brand-teal font-bold">
                                {cov.coverage_breadth_10x !== undefined ? cov.coverage_breadth_10x.toFixed(2) : '0.00'}%
                              </td>
                              <td className="p-4 font-mono text-gray-300">
                                {cons.length !== undefined ? cons.length.toLocaleString() : 0} bp
                              </td>
                              <td className={`p-4 font-mono font-bold ${
                                cons.n_content > 5 ? 'text-brand-rose' : 'text-gray-400'
                              }`}>
                                {cons.n_content !== undefined ? cons.n_content.toFixed(2) : '0.00'}%
                              </td>
                              <td className="p-4 font-mono text-gray-300">
                                {totalSNPs}
                              </td>
                              <td className="p-4 text-center">
                                <div className="inline-flex justify-center">
                                  {status === 'PASS' && (
                                    <span className="flex items-center gap-1 text-[10px] font-bold text-brand-emerald bg-brand-emerald/15 px-2 py-0.5 rounded-full border border-brand-emerald/20">
                                      <CheckCircle2 className="w-3 h-3" /> PASS
                                    </span>
                                  )}
                                  {status === 'WARN' && (
                                    <span className="flex items-center gap-1 text-[10px] font-bold text-brand-amber bg-brand-amber/15 px-2 py-0.5 rounded-full border border-brand-amber/20">
                                      <AlertTriangle className="w-3 h-3" /> WARN
                                    </span>
                                  )}
                                  {status === 'FAIL' && (
                                    <span className="flex items-center gap-1 text-[10px] font-bold text-brand-rose bg-brand-rose/15 px-2 py-0.5 rounded-full border border-brand-rose/20">
                                      <XCircle className="w-3 h-3" /> FAIL
                                    </span>
                                  )}
                                </div>
                              </td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </div>
              </div>
            )}

             {/* TAB CONTENT: SAMPLE DETAILS */}
            {activeTab === 'details' && (() => {
              const refs = getReferencesList();
              const activeCov = (selectedRefName && activeSample.coverage && activeSample.coverage[selectedRefName]) ? activeSample.coverage[selectedRefName] : (activeSample.coverage || {});
              const activeCons = (selectedRefName && activeSample.consensus && activeSample.consensus[selectedRefName]) ? activeSample.consensus[selectedRefName] : (activeSample.consensus || {});
              const activeRefLength = (selectedRefName && refs.find(r => r.name === selectedRefName)) ? refs.find(r => r.name === selectedRefName).length : (reportData.reference?.length || 0);
              const activeGenes = (selectedRefName && reportData.genes) ? (reportData.genes[selectedRefName] || []) : (reportData.genes || []);
              const activeVariants = (selectedRefName && activeSample.variants) ? activeSample.variants.filter(v => v.reference === selectedRefName || !v.reference) : (activeSample.variants || []);

              return (
                <div className="flex flex-col gap-6">
                  
                  {/* Sample Selector row */}
                  <div className="flex flex-col lg:flex-row lg:items-center justify-between gap-4 p-4 glass-panel rounded-2xl">
                    <div className="flex flex-wrap items-center gap-4">
                      <div className="flex items-center gap-3">
                        <span className="text-xs text-gray-400 font-bold uppercase tracking-wider font-outfit">Active Sample:</span>
                        <select
                          value={selectedSampleId}
                          onChange={(e) => {
                            const newSampleId = e.target.value;
                            setSelectedSampleId(newSampleId);
                            const newSample = reportData.samples.find(s => s.sample_id === newSampleId);
                            const bestRef = getBestReference(newSample);
                            if (bestRef) setSelectedRefName(bestRef.name);
                          }}
                          className="bg-dark-bg border border-dark-border text-white text-sm font-normal font-outfit rounded-xl px-4 py-2 focus:outline-none focus:border-brand-teal cursor-pointer"
                        >
                          {[...reportData.samples]
                            .sort((a, b) => a.sample_id.localeCompare(b.sample_id, undefined, { numeric: true, sensitivity: 'base' }))
                            .map(s => (
                              <option key={s.sample_id} value={s.sample_id} className="bg-dark-card font-normal">{s.sample_id}</option>
                            ))
                          }
                        </select>
                      </div>

                      {/* Strain selector - only shown if multi-reference run */}
                      {refs.length > 1 && (
                        <div className="flex items-center gap-3">
                          <span className="text-xs text-gray-400 font-bold uppercase tracking-wider font-outfit">Reference Strain:</span>
                          <select
                            value={selectedRefName}
                            onChange={(e) => setSelectedRefName(e.target.value)}
                            className="bg-dark-bg border border-dark-border text-white text-sm font-normal font-outfit rounded-xl px-4 py-2 focus:outline-none focus:border-brand-teal cursor-pointer"
                          >
                            {refs.map(r => {
                              const cov = activeSample.coverage && activeSample.coverage[r.name];
                              const breadth = cov ? cov.coverage_breadth_10x : 0;
                              return (
                                <option key={r.name} value={r.name} className="bg-dark-card font-normal">
                                  {r.name} ({breadth.toFixed(1)}% cov)
                                </option>
                              );
                            })}
                          </select>
                        </div>
                      )}
                    </div>

                    <div className="flex flex-wrap items-center gap-x-4 gap-y-1 text-xs">
                      <span className="text-gray-400">Mapped Reads:</span>
                      <strong className="text-brand-indigo font-mono text-sm">
                        {(() => {
                          const mr = getMappedReadsForSampleRef(activeSample, selectedRefName);
                          const filtered = activeSample.qc?.trimmed_reads || 0;
                          const pct = filtered > 0 ? ((mr / filtered) * 100).toFixed(1) : '0.0';
                          return `${mr.toLocaleString()} (${pct}% of ${filtered.toLocaleString()} filtered reads)`;
                        })()}
                      </strong>
                      <span className="text-gray-600">|</span>
                      <span className="text-gray-400">Mean Alignment Depth:</span>
                      <strong className="text-brand-teal font-mono text-sm">{activeCov.mean_depth !== undefined ? activeCov.mean_depth : 0}x</strong>
                      <span className="text-gray-600">|</span>
                      <span className="text-gray-400">GC:</span>
                      <strong className="text-white font-mono">{activeCons.gc_content !== undefined ? activeCons.gc_content : 0}%</strong>
                    </div>
                  </div>

                  {activeSample.metadata && Object.keys(activeSample.metadata).length > 0 && (
                    <div className="flex flex-wrap items-center gap-x-6 gap-y-2 text-xs bg-slate-100/80 border border-slate-200 rounded-2xl p-4 shadow-sm">
                      {Object.entries(activeSample.metadata)
                        .filter(([k]) => k !== 'sample_id')
                        .map(([key, val]) => (
                          <span key={key} className="font-outfit">
                            <span className="text-[10px] text-slate-500 uppercase font-bold tracking-wider mr-2">{key.replace('_', ' ')}:</span>
                            <strong className="text-slate-800 font-semibold">{val}</strong>
                          </span>
                        ))
                      }
                    </div>
                  )}

                  {/* 1. Interactive Coverage Chart */}
                  <CoverageChart 
                    coverageData={activeCov.coverage_profile} 
                    minCoverageThreshold={minCovThreshold}
                  />

                  {/* 2. Gene Browser with variants markers */}
                  <GeneMap 
                    genes={activeGenes} 
                    variants={activeVariants} 
                    refLength={activeRefLength} 
                  />

                  {/* 3. Detailed VCF Variant Table */}
                  <VariantTable 
                    variants={activeVariants} 
                    sampleId={activeSample.sample_id} 
                  />
                </div>
              );
            })()}

            {/* TAB CONTENT: MUTATIONAL SHARING */}
            {activeTab === 'matrix' && (
              <MutationMatrix 
                samples={reportData.samples} 
                selectedRefName={selectedRefName} 
                setSelectedRefName={setSelectedRefName}
                refs={getReferencesList()}
              />
            )}

            {/* TAB CONTENT: PHYLOGENETIC TREE */}
            {activeTab === 'tree' && (
              <PhylogeneticTree 
                samples={reportData.samples} 
                selectedRefName={selectedRefName} 
                setSelectedRefName={setSelectedRefName}
                refs={getReferencesList()}
                distanceMatrices={reportData.distance_matrices}
              />
            )}
            {/* TAB CONTENT: GENOMICS PDF REPORT */}
            {activeTab === 'report' && (() => {
              const refs = getReferencesList();
              const pdfRefName = selectedRefName || (refs.length > 0 ? refs[0].name : '');
              const pdfCov = (pdfRefName && activeSample.coverage && activeSample.coverage[pdfRefName]) ? activeSample.coverage[pdfRefName] : (activeSample.coverage || {});
              const pdfCons = (pdfRefName && activeSample.consensus && activeSample.consensus[pdfRefName]) ? activeSample.consensus[pdfRefName] : (activeSample.consensus || {});
              const pdfRefLabel = reportData.reference ? reportData.reference.name : pdfRefName;
              const coDetections = getCoDetections(activeSample);

              return (
                <div className="flex flex-col gap-6">
                  
                  {/* PDF generation options */}
                  <div className="flex items-center justify-between p-4 glass-panel rounded-2xl no-print">
                    <div>
                      <h4 className="text-sm font-semibold font-outfit text-white">Genomic Epidemiological PDF Generator</h4>
                      <p className="text-xs text-gray-400">Review print layout below and click export to generate a clean PDF.</p>
                    </div>
                    
                    <div className="flex items-center gap-3">
                      <select
                        value={selectedSampleId}
                        onChange={(e) => setSelectedSampleId(e.target.value)}
                        className="bg-dark-bg border border-dark-border text-white text-xs font-normal rounded-xl px-3 py-1.5 focus:outline-none cursor-pointer"
                      >
                        {[...reportData.samples]
                          .sort((a, b) => a.sample_id.localeCompare(b.sample_id, undefined, { numeric: true, sensitivity: 'base' }))
                          .map(s => (
                            <option key={s.sample_id} value={s.sample_id} className="bg-dark-card font-normal">{s.sample_id}</option>
                          ))
                        }
                      </select>

                      <button
                        onClick={() => window.print()}
                        className="flex items-center gap-1.5 px-3 py-1.5 bg-brand-teal text-dark-bg text-xs font-bold rounded-xl hover:opacity-90 active:scale-95 transition-all shadow-glow-teal"
                      >
                        <Printer className="w-4 h-4" />
                        Print / Save PDF
                      </button>
                    </div>
                  </div>

                  {/* Report preview container */}
                  <div className="p-8 bg-white text-slate-800 rounded-3xl max-w-4xl mx-auto shadow-2xl flex flex-col gap-8 border border-gray-200 print-card">
                    
                    {/* Report Header */}
                    <div className="flex items-start justify-between border-b border-gray-200 pb-6">
                      <div className="flex flex-col gap-1.5">
                        <span className="text-[10px] font-bold uppercase tracking-wider text-teal-600 font-sans">Genomic Epidemiology Report</span>
                        <h2 className="text-2xl font-bold font-outfit text-slate-900 leading-none">Pathogen genomic sequencing</h2>
                        <span className="text-xs text-slate-500 font-mono mt-1">Report ID: VMR-{activeSample.sample_id}-{Date.now().toString().slice(-6)}</span>
                      </div>
                      <div className="flex flex-col items-end text-right gap-1 text-xs">
                        <span className="font-semibold text-slate-700">Viromics Lab Platform</span>
                        <span className="text-slate-500">Date: {new Date().toLocaleDateString()}</span>
                        <span className="text-slate-500">Run ID: {reportData.run_id}</span>
                      </div>
                    </div>

                    {/* Summary / Patient & sample Info */}
                    <div className="grid grid-cols-2 md:grid-cols-3 gap-6 bg-slate-50 p-5 rounded-2xl border border-slate-100">
                      <div className="flex flex-col gap-1">
                        <span className="text-[10px] uppercase font-bold text-slate-400">Sample Identifier</span>
                        <strong className="text-sm font-semibold text-slate-800">{activeSample.sample_id}</strong>
                      </div>
                      <div className="flex flex-col gap-1">
                        <span className="text-[10px] uppercase font-bold text-slate-400">Reference Genome / Target</span>
                        <strong className="text-sm font-semibold text-slate-800 font-mono">
                          {coDetections.length > 1 ? `Co-detection: ${coDetections.join(' + ')}` : pdfRefLabel}
                        </strong>
                      </div>
                      <div className="flex flex-col gap-1">
                        <span className="text-[10px] uppercase font-bold text-slate-400">Consensus Length</span>
                        <strong className="text-sm font-semibold text-slate-800 font-mono">
                          {coDetections.length > 1 
                            ? `${coDetections.map(c => {
                                const cons = activeSample.consensus && activeSample.consensus[c];
                                return cons ? `${cons.length} bp (${c.split('_')[0]})` : '';
                              }).filter(Boolean).join(' | ')}`
                            : `${pdfCons.length !== undefined ? pdfCons.length.toLocaleString() : 0} bp`
                          }
                        </strong>
                      </div>
                    </div>

                    {/* QC metrics summary */}
                    <div>
                      <h3 className="text-xs font-bold text-slate-900 border-b border-slate-200 pb-2 mb-4 font-outfit">SEQUENCING QUALITY ASSURANCE</h3>
                      {coDetections.length > 1 ? (
                        <div className="overflow-hidden border border-slate-100 rounded-xl bg-slate-50/50">
                          <table className="w-full text-left text-xs border-collapse">
                            <thead>
                              <tr className="bg-slate-100 text-slate-600 font-semibold border-b border-slate-200">
                                <th className="p-3">Detected Strain</th>
                                <th className="p-3">Mapped Reads</th>
                                <th className="p-3">Mean Depth</th>
                                <th className="p-3">Coverage Breadth (≥{minCovThreshold}x)</th>
                                <th className="p-3">Consensus N Rate</th>
                                <th className="p-3 text-center">Status</th>
                              </tr>
                            </thead>
                            <tbody>
                              {coDetections.map(c => {
                                const cov = activeSample.coverage && activeSample.coverage[c] ? activeSample.coverage[c] : {};
                                const cons = activeSample.consensus && activeSample.consensus[c] ? activeSample.consensus[c] : {};
                                const status = getSampleStatusForRef(activeSample, c);
                                const mr = getMappedReadsForSampleRef(activeSample, c);
                                const filtered = activeSample.qc?.trimmed_reads || 0;
                                const pct = filtered > 0 ? ((mr / filtered) * 100).toFixed(1) : '0.0';
                                return (
                                  <tr key={c} className="border-b border-slate-100 text-slate-700">
                                    <td className="p-3 font-semibold text-slate-800">{c}</td>
                                    <td className="p-3 font-mono font-medium text-slate-900">
                                      {mr.toLocaleString()} ({pct}% of {filtered.toLocaleString()})
                                    </td>
                                    <td className="p-3 font-mono">{cov.mean_depth !== undefined ? cov.mean_depth : 0}x</td>
                                    <td className="p-3 font-mono text-teal-600 font-bold">{cov.coverage_breadth_10x !== undefined ? cov.coverage_breadth_10x.toFixed(2) : '0.00'}%</td>
                                    <td className={`p-3 font-mono font-bold ${cons.n_content > 5 ? 'text-red-500' : 'text-slate-500'}`}>
                                      {cons.n_content !== undefined ? cons.n_content.toFixed(2) : '0.00'}%
                                    </td>
                                    <td className="p-3 text-center">
                                      <span className={`inline-block px-2 py-0.5 rounded-full text-[9px] font-bold ${
                                        status === 'PASS' ? 'bg-emerald-100 text-emerald-700' :
                                        status === 'WARN' ? 'bg-amber-100 text-amber-700' : 'bg-red-100 text-red-700'
                                      }`}>
                                        {status}
                                      </span>
                                    </td>
                                  </tr>
                                );
                              })}
                            </tbody>
                          </table>
                        </div>
                      ) : (
                        <div className="grid grid-cols-2 md:grid-cols-5 gap-4 text-center">
                          <div className="p-3 border border-slate-100 rounded-xl bg-slate-50">
                            <span className="text-[9px] uppercase font-bold text-slate-400 block">Total Raw Reads</span>
                            <strong className="text-md text-slate-800 font-semibold font-mono">{activeSample.qc.raw_reads.toLocaleString()}</strong>
                          </div>
                          <div className="p-3 border border-slate-100 rounded-xl bg-slate-50 flex flex-col justify-center items-center">
                            <span className="text-[9px] uppercase font-bold text-slate-400 block">Mapped Reads</span>
                            <strong className="text-md text-brand-indigo font-semibold font-mono leading-tight">
                              {(() => {
                                const mr = getMappedReadsForSampleRef(activeSample, pdfRefName);
                                return mr.toLocaleString();
                              })()}
                            </strong>
                            <span className="text-[9px] text-slate-500 font-mono mt-0.5 leading-none block">
                              {(() => {
                                const mr = getMappedReadsForSampleRef(activeSample, pdfRefName);
                                const filtered = activeSample.qc?.trimmed_reads || 0;
                                const pct = filtered > 0 ? ((mr / filtered) * 100).toFixed(1) : '0.0';
                                return `${pct}% of filtered`;
                              })()}
                            </span>
                          </div>
                          <div className="p-3 border border-slate-100 rounded-xl bg-slate-50">
                            <span className="text-[9px] uppercase font-bold text-slate-400 block">Mean Alignment Depth</span>
                            <strong className="text-md text-teal-600 font-semibold font-mono">{pdfCov.mean_depth !== undefined ? pdfCov.mean_depth : 0}x</strong>
                          </div>
                          <div className="p-3 border border-slate-100 rounded-xl bg-slate-50">
                            <span className="text-[9px] uppercase font-bold text-slate-400 block">Coverage Breadth (≥{minCovThreshold}x)</span>
                            <strong className="text-md text-teal-600 font-semibold font-mono">{pdfCov.coverage_breadth_10x !== undefined ? pdfCov.coverage_breadth_10x : 0}%</strong>
                          </div>
                          <div className="p-3 border border-slate-100 rounded-xl bg-slate-50">
                            <span className="text-[9px] uppercase font-bold text-slate-400 block">Consensus N Rate</span>
                            <strong className={`text-md font-semibold font-mono ${
                              pdfCons.n_content > 5 ? 'text-red-500' : 'text-slate-800'
                            }`}>{pdfCons.n_content !== undefined ? pdfCons.n_content : 0}%</strong>
                          </div>
                        </div>
                      )}
                    </div>

                    {/* Key Variants Called */}
                    <div>
                      <h3 className="text-xs font-bold text-slate-900 border-b border-slate-200 pb-2 mb-4 font-outfit">IDENTIFIED MUTATIONS</h3>
                      <div className="overflow-hidden border border-slate-200 rounded-xl">
                        <table className="w-full text-left text-xs border-collapse">
                          <thead>
                            <tr className="bg-slate-50 border-b border-slate-200 text-slate-500">
                              {refs.length > 1 && <th className="p-3">Reference Strain</th>}
                              <th className="p-3">Position</th>
                              <th className="p-3">Nucleotide Swap</th>
                              <th className="p-3">Gene Region</th>
                              <th className="p-3">Protein Change</th>
                              <th className="p-3">Type</th>
                              <th className="p-3 text-right">Allele Frequency</th>
                            </tr>
                          </thead>
                          <tbody>
                            {activeSample.variants && activeSample.variants.length > 0 ? (
                              activeSample.variants.map((v, i) => (
                                <tr key={i} className="border-b border-slate-100 text-slate-700">
                                  {refs.length > 1 && (
                                    <td className="p-3 font-semibold text-slate-800 font-outfit">
                                      {v.reference ? v.reference.split('_')[0] : 'All'}
                                    </td>
                                  )}
                                  <td className="p-3 font-mono font-medium">{v.position}</td>
                                  <td className="p-3 font-mono font-bold text-slate-900">{v.ref} &gt; {v.alt}</td>
                                  <td className="p-3">
                                    <span className="px-2 py-0.5 rounded text-[10px] font-semibold bg-indigo-50 text-indigo-700 border border-indigo-100">
                                      {v.gene}
                                    </span>
                                  </td>
                                  <td className="p-3 font-medium">
                                    {v.aa_change ? (
                                      v.aa_change.includes('(syn)') ? (
                                        <span className="text-slate-400 font-normal">{v.aa_change}</span>
                                      ) : (
                                        <span className="font-bold text-red-600 font-mono">p.{v.aa_change}</span>
                                      )
                                    ) : (
                                      <span className="text-slate-400 font-normal">-</span>
                                    )}
                                  </td>
                                  <td className="p-3">{v.type}</td>
                                  <td className="p-3 font-mono text-right text-teal-600 font-bold">{(v.frequency * 100).toFixed(1)}%</td>
                                </tr>
                              ))
                            ) : (
                              <tr>
                                <td colSpan={refs.length > 1 ? 7 : 6} className="p-4 text-center text-slate-400">No variants identified in this run.</td>
                              </tr>
                            )}
                          </tbody>
                        </table>
                      </div>
                    </div>

                    {/* Signatures */}
                    <div className="flex justify-between items-center border-t border-slate-200 pt-8 mt-4 text-xs">
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

          </div>
        )}

      </main>
    </div>
  );
}
