import React, { useState, useEffect } from 'react';
import { Search, Filter, Download, ArrowUpDown, ChevronLeft, ChevronRight } from 'lucide-react';
import { formatConsequences, formatMetric } from '../utils/reportSchema';
import { formatVariantsCsv, triggerBrowserDownload } from '../utils/exporters';

export default function VariantTable({ variants, sampleId }) {
  const [searchTerm, setSearchTerm] = useState('');
  const [geneFilter, setGeneFilter] = useState('ALL');
  const [typeFilter, setTypeFilter] = useState('ALL');
  const [sortField, setSortField] = useState('position');
  const [sortAsc, setSortAsc] = useState(true);

  const [currentPage, setCurrentPage] = useState(1);
  const [itemsPerPage, setItemsPerPage] = useState(10);

  useEffect(() => {
    setCurrentPage(1);
  }, [sampleId, searchTerm, geneFilter, typeFilter]);

  if (!variants || variants.length === 0) {
    return (
      <div className="p-8 text-center text-slate-500 glass-panel rounded-2xl">
        <p>No variants identified for sample <strong className="text-slate-900">{sampleId}</strong>.</p>
      </div>
    );
  }

  const uniqueGenes = [
    'ALL',
    ...new Set(
      variants.flatMap((v) => {
        const cons = formatConsequences(v);
        return cons.length > 0 ? cons.map((c) => c.gene) : [v.gene];
      }).filter(Boolean)
    ),
  ];

  const handleSort = (field) => {
    if (sortField === field) {
      setSortAsc(!sortAsc);
    } else {
      setSortField(field);
      setSortAsc(true);
    }
  };

  const filteredVariants = variants.filter((v) => {
    const cons = formatConsequences(v);
    const geneStr = v.gene || '';
    const aaStr = v.aa_change || '';
    const statusStr = v.status || '';
    const q = searchTerm.toLowerCase();
    const matchesSearch =
      !q ||
      v.position.toString().includes(searchTerm) ||
      geneStr.toLowerCase().includes(q) ||
      aaStr.toLowerCase().includes(q) ||
      statusStr.toLowerCase().includes(q) ||
      cons.some((c) => `${c.gene} ${c.label}`.toLowerCase().includes(q)) ||
      `${v.ref}>${v.alt}`.toLowerCase().includes(q);

    const matchesGene =
      geneFilter === 'ALL' ||
      v.gene === geneFilter ||
      cons.some((c) => c.gene === geneFilter);
    const matchesType = typeFilter === 'ALL' || v.type === typeFilter;

    return matchesSearch && matchesGene && matchesType;
  });

  const sortedVariants = [...filteredVariants].sort((a, b) => {
    let valA = a[sortField];
    let valB = b[sortField];

    if (valA === null || valA === undefined) return 1;
    if (valB === null || valB === undefined) return -1;

    if (typeof valA === 'string') {
      valA = valA.toLowerCase();
      valB = String(valB).toLowerCase();
    }

    if (valA < valB) return sortAsc ? -1 : 1;
    if (valA > valB) return sortAsc ? 1 : -1;
    return 0;
  });

  const currentItemsPerPage = itemsPerPage === 'ALL' ? sortedVariants.length : Number(itemsPerPage);
  const totalPages = Math.ceil(sortedVariants.length / currentItemsPerPage) || 1;
  const startIndex = (currentPage - 1) * currentItemsPerPage;
  const paginatedVariants = sortedVariants.slice(startIndex, startIndex + currentItemsPerPage);

  const exportToCSV = () => {
    const csv = formatVariantsCsv(sortedVariants, sampleId);
    triggerBrowserDownload(csv, `viromics_${sampleId}_variants.csv`, 'text/csv;charset=utf-8;');
  };

  return (
    <div className="p-6 glass-panel rounded-2xl flex flex-col gap-6">
      <div className="flex flex-col xl:flex-row xl:items-center justify-between gap-4">
        <div>
          <h3 className="text-lg font-semibold text-slate-900">Variant Call Details</h3>
          <p className="text-xs text-slate-500">
            High-confidence SNPs &amp; indels computed from FORMAT/AD &amp; FORMAT/DP with consensus/subclonal status
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-3">
          <div className="relative">
            <Search className="w-4 h-4 absolute left-3 top-2.5 text-slate-400" />
            <input
              type="text"
              placeholder="Search pos, gene, AA..."
              value={searchTerm}
              onChange={(e) => {
                setSearchTerm(e.target.value);
                setCurrentPage(1);
              }}
              className="pl-9 pr-4 py-2 w-48 text-xs rounded-xl bg-slate-50 border border-slate-200 text-slate-900 focus:outline-none focus:border-red-500 transition-all"
            />
          </div>

          <div className="flex items-center gap-1 bg-slate-50 border border-slate-200 rounded-xl px-2.5 py-1.5">
            <Filter className="w-3.5 h-3.5 text-red-600" />
            <select
              value={geneFilter}
              onChange={(e) => {
                setGeneFilter(e.target.value);
                setCurrentPage(1);
              }}
              className="text-xs bg-transparent text-slate-700 focus:outline-none border-none cursor-pointer"
            >
              {uniqueGenes.map((g) => (
                <option key={g} value={g}>
                  {g === 'ALL' ? 'All Genes' : g}
                </option>
              ))}
            </select>
          </div>

          <select
            value={typeFilter}
            onChange={(e) => {
              setTypeFilter(e.target.value);
              setCurrentPage(1);
            }}
            className="text-xs bg-slate-50 border border-slate-200 text-slate-700 rounded-xl px-3 py-1.5 focus:outline-none cursor-pointer"
          >
            <option value="ALL">All Types</option>
            <option value="Coding">Coding Only</option>
            <option value="Non-coding">Non-Coding Only</option>
          </select>

          <button
            onClick={exportToCSV}
            className="flex items-center gap-1.5 px-3 py-1.5 bg-red-600 text-white text-xs font-semibold rounded-xl hover:bg-red-700 active:scale-95 transition-all"
          >
            <Download className="w-3.5 h-3.5" />
            Export CSV
          </button>
        </div>
      </div>

      <div className="overflow-x-auto border border-slate-200 rounded-xl bg-white max-w-full">
        <table className="w-full text-left text-xs border-collapse">
          <thead>
            <tr className="border-b border-slate-200 text-slate-700 bg-slate-100">
              <th className="p-3.5 cursor-pointer hover:text-slate-900" onClick={() => handleSort('position')}>
                <div className="flex items-center gap-1.5">
                  Position <ArrowUpDown className="w-3 h-3 text-red-600" />
                </div>
              </th>
              <th className="p-3.5">Allele &amp; Consequence</th>
              <th className="p-3.5 cursor-pointer hover:text-slate-900" onClick={() => handleSort('gene')}>
                <div className="flex items-center gap-1.5">
                  Gene Region <ArrowUpDown className="w-3 h-3 text-red-600" />
                </div>
              </th>
              <th className="p-3.5">Class / Consensus</th>
              <th className="p-3.5 cursor-pointer hover:text-slate-900 text-right" onClick={() => handleSort('depth')}>
                <div className="flex items-center justify-end gap-1.5">
                  Depth <ArrowUpDown className="w-3 h-3 text-red-600" />
                </div>
              </th>
              <th className="p-3.5 cursor-pointer hover:text-slate-900 text-right" onClick={() => handleSort('frequency')}>
                <div className="flex items-center justify-end gap-1.5">
                  Allele Freq <ArrowUpDown className="w-3 h-3 text-red-600" />
                </div>
              </th>
              <th className="p-3.5 cursor-pointer hover:text-slate-900 text-right" onClick={() => handleSort('quality')}>
                <div className="flex items-center justify-end gap-1.5">
                  Quality <ArrowUpDown className="w-3 h-3 text-red-600" />
                </div>
              </th>
            </tr>
          </thead>
          <tbody>
            {paginatedVariants.length > 0 ? (
              paginatedVariants.map((v, idx) => {
                const cons = formatConsequences(v);
                return (
                  <tr key={idx} className="border-b border-slate-100 hover:bg-slate-50 transition-colors">
                    <td className="p-3.5 font-mono font-medium text-slate-900">
                      {v.position.toLocaleString()}
                    </td>
                    <td className="p-3.5 font-mono">
                      <div className="flex flex-col gap-1">
                        <div>
                          <span className="text-slate-500">{v.ref}</span>
                          <span className="text-indigo-600 font-bold mx-1.5">&gt;</span>
                          <span className="text-red-600 font-bold">{v.alt}</span>
                        </div>
                        {cons.length > 0 && (
                          <div className="flex flex-wrap gap-1">
                            {cons.map((c, cIdx) => (
                              <span
                                key={cIdx}
                                title={c.product || c.status || c.effect}
                                className={`text-[10px] px-1.5 py-0.5 rounded ${
                                  c.kind === 'unavailable'
                                    ? 'bg-slate-100 text-slate-500 italic'
                                    : c.kind === 'synonymous'
                                    ? 'bg-slate-100 text-slate-500'
                                    : c.kind === 'frameshift' || c.kind === 'indel'
                                    ? 'bg-amber-50 text-amber-700 font-semibold'
                                    : 'bg-red-50 text-red-700 font-bold'
                                }`}
                              >
                                {cons.length > 1 ? `${c.gene}: ` : ''}
                                {c.label}
                              </span>
                            ))}
                          </div>
                        )}
                      </div>
                    </td>
                    <td className="p-3.5">
                      <span className="px-2 py-0.5 rounded text-[10px] font-bold bg-indigo-50 text-indigo-700 border border-indigo-200">
                        {v.gene}
                      </span>
                    </td>
                    <td className="p-3.5">
                      <div className="flex flex-wrap gap-1.5 items-center">
                        <span className={`px-2 py-0.5 rounded-full text-[10px] border font-medium ${
                          v.type === 'Coding'
                            ? 'bg-red-50 text-red-700 border-red-200'
                            : 'bg-slate-100 text-slate-600 border-slate-200'
                        }`}>
                          {v.type}
                        </span>
                        {v.is_consensus !== undefined && (
                          <span className={`px-2 py-0.5 rounded-full text-[10px] font-semibold ${
                            v.is_consensus
                              ? 'bg-emerald-50 text-emerald-700 border border-emerald-200'
                              : 'bg-amber-50 text-amber-700 border border-amber-200'
                          }`}>
                            {v.is_consensus ? 'Consensus' : 'Subclonal'}
                          </span>
                        )}
                      </div>
                    </td>
                    <td className="p-3.5 font-mono text-right text-slate-700">{v.depth}x</td>
                    <td className="p-3.5 font-mono text-right text-red-600 font-semibold">
                      {(v.frequency * 100).toFixed(1)}%
                    </td>
                    <td className="p-3.5 font-mono text-right text-slate-700">
                      {formatMetric(v.quality, { decimals: 1 })}
                    </td>
                  </tr>
                );
              })
            ) : (
              <tr>
                <td colSpan="7" className="p-8 text-center text-slate-500 font-medium">
                  No variants found matching the filters.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {sortedVariants.length > 0 && (
        <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4 border-t border-slate-200 pt-4 text-xs text-slate-500">
          <div className="flex flex-wrap items-center gap-4">
            <span>
              Showing <strong className="text-slate-900">{startIndex + 1}</strong> to{' '}
              <strong className="text-slate-900">{Math.min(startIndex + currentItemsPerPage, sortedVariants.length)}</strong> of{' '}
              <strong className="text-slate-900">{sortedVariants.length}</strong> variants
            </span>

            <div className="flex items-center gap-1.5">
              <span>Rows:</span>
              <select
                value={itemsPerPage}
                onChange={(e) => {
                  const val = e.target.value;
                  setItemsPerPage(val === 'ALL' ? 'ALL' : Number(val));
                  setCurrentPage(1);
                }}
                className="bg-slate-50 border border-slate-200 text-slate-800 text-[11px] rounded-lg px-2 py-0.5 focus:outline-none cursor-pointer"
              >
                <option value={10}>10</option>
                <option value={50}>50</option>
                <option value={200}>200</option>
                <option value="ALL">All</option>
              </select>
            </div>
          </div>

          {totalPages > 1 && (
            <div className="flex items-center gap-2">
              <button
                onClick={() => setCurrentPage((prev) => Math.max(1, prev - 1))}
                disabled={currentPage === 1}
                className="p-1.5 rounded-lg border border-slate-200 bg-white text-slate-700 hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed"
              >
                <ChevronLeft className="w-4 h-4" />
              </button>
              <span className="px-2">
                Page <strong className="text-slate-900">{currentPage}</strong> of{' '}
                <strong className="text-slate-900">{totalPages}</strong>
              </span>
              <button
                onClick={() => setCurrentPage((prev) => Math.min(totalPages, prev + 1))}
                disabled={currentPage === totalPages}
                className="p-1.5 rounded-lg border border-slate-200 bg-white text-slate-700 hover:bg-slate-50 disabled:opacity-40 disabled:cursor-not-allowed"
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
