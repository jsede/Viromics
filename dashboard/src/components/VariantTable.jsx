import React, { useState } from 'react';
import { Search, Filter, Download, ArrowUpDown, ChevronLeft, ChevronRight } from 'lucide-react';

export default function VariantTable({ variants, sampleId }) {
  const [searchTerm, setSearchTerm] = useState('');
  const [geneFilter, setGeneFilter] = useState('ALL');
  const [typeFilter, setTypeFilter] = useState('ALL');
  const [sortField, setSortField] = useState('position');
  const [sortAsc, setSortAsc] = useState(true);
  
  // Pagination
  const [currentPage, setCurrentPage] = useState(1);
  const [itemsPerPage, setItemsPerPage] = useState(10);

  if (!variants || variants.length === 0) {
    return (
      <div className="p-8 text-center text-gray-400 glass-panel rounded-2xl">
        <p>No variants identified for sample <strong className="text-white">{sampleId}</strong>.</p>
      </div>
    );
  }

  // Get list of unique genes for filtering dropdown
  const uniqueGenes = ['ALL', ...new Set(variants.map(v => v.gene))];

  // Sorting handler
  const handleSort = (field) => {
    if (sortField === field) {
      setSortAsc(!sortAsc);
    } else {
      setSortField(field);
      setSortAsc(true);
    }
  };

  // Filter & Search Logic
  const filteredVariants = variants.filter(v => {
    const matchesSearch = 
      v.position.toString().includes(searchTerm) ||
      v.gene.toLowerCase().includes(searchTerm.toLowerCase()) ||
      `${v.ref}>${v.alt}`.toLowerCase().includes(searchTerm.toLowerCase());
      
    const matchesGene = geneFilter === 'ALL' || v.gene === geneFilter;
    const matchesType = typeFilter === 'ALL' || v.type === typeFilter;
    
    return matchesSearch && matchesGene && matchesType;
  });

  // Sort Logic
  const sortedVariants = [...filteredVariants].sort((a, b) => {
    let valA = a[sortField];
    let valB = b[sortField];
    
    if (typeof valA === 'string') {
      valA = valA.toLowerCase();
      valB = valB.toLowerCase();
    }
    
    if (valA < valB) return sortAsc ? -1 : 1;
    if (valA > valB) return sortAsc ? 1 : -1;
    return 0;
  });

  // Paginated slices
  const currentItemsPerPage = itemsPerPage === 'ALL' ? sortedVariants.length : Number(itemsPerPage);
  const totalPages = Math.ceil(sortedVariants.length / currentItemsPerPage) || 1;
  const startIndex = (currentPage - 1) * currentItemsPerPage;
  const paginatedVariants = sortedVariants.slice(startIndex, startIndex + currentItemsPerPage);

  // CSV Export utility
  const exportToCSV = () => {
    const headers = ["Position", "Ref", "Alt", "Gene", "Type", "Depth", "Frequency", "Quality"];
    const csvRows = [headers.join(",")];
    
    sortedVariants.forEach(v => {
      csvRows.push([
        v.position,
        v.ref,
        v.alt,
        `"${v.gene}"`,
        v.type,
        v.depth,
        v.frequency,
        v.quality
      ].join(","));
    });
    
    const blob = new Blob([csvRows.join("\n")], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement("a");
    link.setAttribute("href", url);
    link.setAttribute("download", `viromics_${sampleId}_variants.csv`);
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  };

  const getBadgeStyle = (type) => {
    if (type === 'Coding') return 'bg-brand-teal/15 text-brand-teal border-brand-teal/20';
    return 'bg-gray-700/35 text-gray-400 border-dark-border';
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      {/* Table Header and Controls */}
      <div className="flex flex-col xl:flex-row xl:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold font-outfit text-white">Variant Call Details</h3>
          <p className="text-xs text-gray-400">
            Browse and filter individual high-confidence SNP and Indel classifications
          </p>
        </div>

        {/* Filters Grid */}
        <div className="flex flex-wrap items-center gap-3">
          {/* Search box */}
          <div className="relative">
            <Search className="w-4 h-4 absolute left-3 top-2.5 text-gray-400" />
            <input
              type="text"
              placeholder="Search pos/gene..."
              value={searchTerm}
              onChange={(e) => { setSearchTerm(e.target.value); setCurrentPage(1); }}
              className="pl-9 pr-4 py-2 w-48 text-xs rounded-xl bg-dark-bg border border-dark-border text-white focus:outline-none focus:border-brand-teal transition-all"
            />
          </div>

          {/* Gene Select */}
          <div className="flex items-center gap-1 bg-dark-bg border border-dark-border rounded-xl px-2.5 py-1.5">
            <Filter className="w-3.5 h-3.5 text-brand-teal" />
            <select
              value={geneFilter}
              onChange={(e) => { setGeneFilter(e.target.value); setCurrentPage(1); }}
              className="text-xs bg-transparent text-gray-300 focus:outline-none border-none cursor-pointer"
            >
              {uniqueGenes.map(g => (
                <option key={g} value={g} className="bg-dark-card">{g === 'ALL' ? 'All Genes' : g}</option>
              ))}
            </select>
          </div>

          {/* Coding / Non-coding Select */}
          <select
            value={typeFilter}
            onChange={(e) => { setTypeFilter(e.target.value); setCurrentPage(1); }}
            className="text-xs bg-dark-bg border border-dark-border text-gray-300 rounded-xl px-3 py-1.5 focus:outline-none cursor-pointer"
          >
            <option value="ALL" className="bg-dark-card">All Types</option>
            <option value="Coding" className="bg-dark-card">Coding Only</option>
            <option value="Non-coding" className="bg-dark-card">Non-Coding Only</option>
          </select>

          {/* Export Button */}
          <button
            onClick={exportToCSV}
            className="flex items-center gap-1.5 px-3 py-1.5 bg-brand-teal text-dark-bg text-xs font-semibold rounded-xl hover:opacity-90 active:scale-95 transition-all shadow-glow-teal"
          >
            <Download className="w-3.5 h-3.5" />
            Export CSV
          </button>
        </div>
      </div>

      {/* Grid Table */}
      <div className="overflow-x-auto border border-dark-border rounded-xl bg-dark-bg/40 max-w-full">
        <table className="w-full text-left text-xs border-collapse">
          <thead>
            <tr className="border-b border-dark-border text-gray-400 bg-dark-bg/60">
              <th 
                className="p-4 cursor-pointer hover:text-white transition-colors"
                onClick={() => handleSort('position')}
              >
                <div className="flex items-center gap-1.5">
                  Position <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th className="p-4">Ref &gt; Alt</th>
              <th 
                className="p-4 cursor-pointer hover:text-white transition-colors"
                onClick={() => handleSort('gene')}
              >
                <div className="flex items-center gap-1.5">
                  Gene Region <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th className="p-4">Classification</th>
              <th 
                className="p-4 cursor-pointer hover:text-white transition-colors text-right"
                onClick={() => handleSort('depth')}
              >
                <div className="flex items-center justify-end gap-1.5">
                  Depth <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th 
                className="p-4 cursor-pointer hover:text-white transition-colors text-right"
                onClick={() => handleSort('frequency')}
              >
                <div className="flex items-center justify-end gap-1.5">
                  Allele Freq <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
              <th 
                className="p-4 cursor-pointer hover:text-white transition-colors text-right"
                onClick={() => handleSort('quality')}
              >
                <div className="flex items-center justify-end gap-1.5">
                  Quality <ArrowUpDown className="w-3 h-3" />
                </div>
              </th>
            </tr>
          </thead>
          <tbody>
            {paginatedVariants.length > 0 ? (
              paginatedVariants.map((v, idx) => (
                <tr 
                  key={idx} 
                  className="border-b border-dark-border/40 hover:bg-dark-hover/40 transition-colors"
                >
                  <td className="p-4 font-mono font-medium text-gray-200">
                    {v.position.toLocaleString()}
                  </td>
                  <td className="p-4 font-mono">
                    <div className="flex flex-col">
                      <div>
                        <span className="text-gray-400">{v.ref}</span>
                        <span className="text-brand-indigo font-bold mx-1.5">&gt;</span>
                        <span className="text-brand-rose font-bold">{v.alt}</span>
                      </div>
                      {v.aa_change && (
                        <span className={`text-[10px] mt-0.5 font-bold ${
                          v.aa_change.includes('(syn)') 
                            ? 'text-gray-500 font-normal' 
                            : 'text-brand-teal'
                        }`}>
                          {v.aa_change.includes('(syn)') ? `${v.aa_change.replace(' (syn)', '')} (syn)` : `p.${v.aa_change}`}
                        </span>
                      )}
                    </div>
                  </td>
                  <td className="p-4">
                    <span className="px-2 py-0.5 rounded text-[10px] font-bold bg-brand-indigo/10 text-brand-indigo border border-brand-indigo/25">
                      {v.gene}
                    </span>
                  </td>
                  <td className="p-4">
                    <span className={`px-2 py-0.5 rounded-full text-[10px] border font-medium ${getBadgeStyle(v.type)}`}>
                      {v.type}
                    </span>
                  </td>
                  <td className="p-4 font-mono text-right text-gray-300">
                    {v.depth}x
                  </td>
                  <td className="p-4 font-mono text-right text-brand-teal font-semibold">
                    {(v.frequency * 100).toFixed(1)}%
                  </td>
                  <td className="p-4 font-mono text-right text-gray-300">
                    {v.quality}
                  </td>
                </tr>
              ))
            ) : (
              <tr>
                <td colSpan="7" className="p-8 text-center text-gray-500 font-medium">
                  No variants found matching the filters.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {/* Pagination controls */}
      {sortedVariants.length > 0 && (
        <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4 border-t border-dark-border/40 pt-4 text-xs text-gray-400">
          <div className="flex flex-wrap items-center gap-4">
            <span>
              Showing <strong className="text-white">{startIndex + 1}</strong> to{' '}
              <strong className="text-white">{Math.min(startIndex + currentItemsPerPage, sortedVariants.length)}</strong> of{' '}
              <strong className="text-white">{sortedVariants.length}</strong> variants
            </span>
            
            {/* Rows Per Page Selector */}
            <div className="flex items-center gap-1.5">
              <span className="text-gray-500">Rows:</span>
              <select
                value={itemsPerPage}
                onChange={(e) => {
                  const val = e.target.value;
                  setItemsPerPage(val === 'ALL' ? 'ALL' : Number(val));
                  setCurrentPage(1);
                }}
                className="bg-dark-bg border border-dark-border text-white text-[10px] font-normal rounded-lg px-2 py-0.5 focus:outline-none cursor-pointer"
              >
                <option value={10} className="bg-dark-card font-normal">10</option>
                <option value={50} className="bg-dark-card font-normal">50</option>
                <option value={200} className="bg-dark-card font-normal">200</option>
                <option value="ALL" className="bg-dark-card font-normal">All</option>
              </select>
            </div>
          </div>

          {totalPages > 1 && (
            <div className="flex items-center gap-2">
              <button
                onClick={() => setCurrentPage(prev => Math.max(1, prev - 1))}
                disabled={currentPage === 1}
                className="p-1.5 rounded-lg border border-dark-border bg-dark-bg hover:bg-dark-hover disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
              >
                <ChevronLeft className="w-4 h-4" />
              </button>
              <span className="px-2">
                Page <strong className="text-white">{currentPage}</strong> of{' '}
                <strong className="text-white">{totalPages}</strong>
              </span>
              <button
                onClick={() => setCurrentPage(prev => Math.min(totalPages, prev + 1))}
                disabled={currentPage === totalPages}
                className="p-1.5 rounded-lg border border-dark-border bg-dark-bg hover:bg-dark-hover disabled:opacity-40 disabled:cursor-not-allowed transition-colors"
              >
                <ChevronRight className="w-4 h-4" />
              </button>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
