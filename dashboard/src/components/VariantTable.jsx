import React, { useState, useEffect } from 'react';
import { Search, Filter, Download, ArrowUpDown, ChevronLeft, ChevronRight } from 'lucide-react';
import { deduplicateConsequences, formatConsequences, formatDepth, formatMetric } from '../utils/reportSchema';
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
      <div className="p-8 text-center text-[#6B7280] glass-panel text-[13px]">
        <p>
          No variants identified for sample <strong className="text-[#111827]">{sampleId}</strong>.
        </p>
      </div>
    );
  }

  const uniqueGenes = [
    'ALL',
    ...new Set(
      variants
        .flatMap((v) => {
          const cons = formatConsequences(v);
          return cons.length > 0 ? cons.map((c) => c.gene) : [v.gene];
        })
        .filter(Boolean)
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

    const matchesGene = geneFilter === 'ALL' || v.gene === geneFilter || cons.some((c) => c.gene === geneFilter);
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
    <div className="p-5 glass-panel flex flex-col gap-5">
      <div className="flex flex-col xl:flex-row xl:items-center justify-between gap-3">
        <div>
          <h3 className="text-[16px] font-semibold text-[#111827]">Variant Call Details</h3>
          <p className="text-[12px] text-[#6B7280] mt-0.5">
            High-confidence SNPs &amp; indels computed from FORMAT/AD &amp; FORMAT/DP with consensus/subclonal status
          </p>
        </div>

        <div className="flex flex-wrap items-center gap-2.5">
          <div className="relative">
            <Search className="w-3.5 h-3.5 absolute left-2.5 top-2.5 text-[#6B7280]" />
            <input
              type="text"
              placeholder="Search pos, gene, AA..."
              value={searchTerm}
              onChange={(e) => {
                setSearchTerm(e.target.value);
                setCurrentPage(1);
              }}
              className="pl-8 pr-3 py-1.5 w-48 text-xs rounded-[4px] bg-white border border-[#CBD2D9] text-[#111827] placeholder-[#6B7280] focus:outline-none focus:border-[#0F5C6E]"
            />
          </div>

          <div className="flex items-center gap-1.5 bg-white border border-[#CBD2D9] rounded-[4px] px-2.5 py-1.5">
            <Filter className="w-3.5 h-3.5 text-[#0F5C6E]" />
            <select
              value={geneFilter}
              onChange={(e) => {
                setGeneFilter(e.target.value);
                setCurrentPage(1);
              }}
              className="text-xs bg-transparent text-[#111827] focus:outline-none border-none cursor-pointer"
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
            className="text-xs bg-white border border-[#CBD2D9] text-[#111827] rounded-[4px] px-2.5 py-1.5 focus:outline-none focus:border-[#0F5C6E] cursor-pointer"
          >
            <option value="ALL">All Types</option>
            <option value="Coding">Coding Only</option>
            <option value="Non-coding">Non-Coding Only</option>
          </select>

          <button
            onClick={exportToCSV}
            className="flex items-center gap-1.5 px-3 py-1.5 bg-[#0F5C6E] hover:bg-[#0B4654] text-white text-xs font-semibold rounded-[4px] transition-colors whitespace-nowrap"
          >
            <Download className="w-3.5 h-3.5" />
            Export CSV
          </button>
        </div>
      </div>

      <div className="overflow-x-auto border border-[#E2E6EA] rounded-[6px] bg-white max-w-full">
        <table className="w-full text-left text-[13px] border-collapse">
          <thead>
            <tr className="border-b border-[#E2E6EA] text-[#4B5563] bg-[#F1F3F5] text-[11px] uppercase tracking-[0.04em] select-none">
              <th className="px-3.5 py-2.5 cursor-pointer hover:text-[#111827] whitespace-nowrap" onClick={() => handleSort('position')}>
                <div className="flex items-center gap-1">
                  Position
                  <ArrowUpDown className={`w-3 h-3 ${sortField === 'position' ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'}`} />
                </div>
              </th>
              <th className="px-3.5 py-2.5">Allele &amp; Consequence</th>
              <th className="px-3.5 py-2.5 cursor-pointer hover:text-[#111827] whitespace-nowrap" onClick={() => handleSort('gene')}>
                <div className="flex items-center gap-1">
                  Gene Region
                  <ArrowUpDown className={`w-3 h-3 ${sortField === 'gene' ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'}`} />
                </div>
              </th>
              <th className="px-3.5 py-2.5 whitespace-nowrap">Class / Consensus</th>
              <th className="px-3.5 py-2.5 cursor-pointer hover:text-[#111827] text-right whitespace-nowrap" onClick={() => handleSort('depth')}>
                <div className="flex items-center justify-end gap-1">
                  Depth
                  <ArrowUpDown className={`w-3 h-3 ${sortField === 'depth' ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'}`} />
                </div>
              </th>
              <th className="px-3.5 py-2.5 cursor-pointer hover:text-[#111827] text-right whitespace-nowrap" onClick={() => handleSort('frequency')}>
                <div className="flex items-center justify-end gap-1">
                  Allele Freq
                  <ArrowUpDown className={`w-3 h-3 ${sortField === 'frequency' ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'}`} />
                </div>
              </th>
              <th className="px-3.5 py-2.5 cursor-pointer hover:text-[#111827] text-right whitespace-nowrap" onClick={() => handleSort('quality')}>
                <div className="flex items-center justify-end gap-1">
                  Quality
                  <ArrowUpDown className={`w-3 h-3 ${sortField === 'quality' ? 'text-[#0F5C6E]' : 'text-[#9CA3AF]'}`} />
                </div>
              </th>
            </tr>
          </thead>
          <tbody>
            {paginatedVariants.length > 0 ? (
              paginatedVariants.map((v, idx) => {
                const dedupCons = deduplicateConsequences(formatConsequences(v));
                return (
                  <tr key={idx} className="border-b border-[#E2E6EA] last:border-none hover:bg-[#F8F9FA] transition-colors">
                    <td className="px-3.5 py-2.5 font-mono font-medium text-[#111827] whitespace-nowrap">
                      {v.position.toLocaleString()}
                    </td>
                    <td className="px-3.5 py-2.5 font-mono">
                      <div className="flex flex-col gap-1">
                        <div className="whitespace-nowrap">
                          <span className="text-[#6B7280]">{v.ref}</span>
                          <span className="text-[#9CA3AF] font-semibold mx-1">&gt;</span>
                          <span className="text-[#0F5C6E] font-bold">{v.alt}</span>
                        </div>
                        {dedupCons.length > 0 && (
                          <div className="flex flex-wrap gap-1">
                            {dedupCons.map((c, cIdx) => (
                              <span
                                key={cIdx}
                                title={c.tooltip || c.label}
                                className={`text-[11px] px-1.5 py-0.5 rounded-[4px] border ${
                                  c.kind === 'unavailable'
                                    ? 'bg-[#F1F3F5] text-[#6B7280] border-[#E2E6EA] italic'
                                    : c.kind === 'synonymous'
                                    ? 'bg-[#F1F3F5] text-[#5C6773] border-[#E2E6EA]'
                                    : c.kind === 'frameshift' || c.kind === 'indel'
                                    ? 'bg-[#FDF6E7] text-[#9A5B06] border-[#EBD29B] font-semibold'
                                    : 'bg-[#FCEEEB] text-[#B43B2E] border-[#EBBDB8] font-semibold'
                                }`}
                              >
                                {dedupCons.length > 1 ? `${c.gene}: ` : ''}
                                {c.label}
                              </span>
                            ))}
                          </div>
                        )}
                      </div>
                    </td>
                    <td className="px-3.5 py-2.5 whitespace-nowrap">
                      <span className="px-2 py-0.5 rounded-[4px] text-[11px] font-semibold bg-[#F1F3F5] text-[#111827] border border-[#E2E6EA] font-mono">
                        {v.gene}
                      </span>
                    </td>
                    <td className="px-3.5 py-2.5 whitespace-nowrap">
                      <div className="flex items-center gap-1.5">
                        <span className="px-2 py-0.5 rounded-[4px] text-[11px] border font-medium bg-[#F8F9FA] text-[#4B5563] border-[#E2E6EA]">
                          {v.type}
                        </span>
                        {v.is_consensus !== undefined && (
                          <span
                            className={`px-2 py-0.5 rounded-[4px] text-[11px] font-semibold ${
                              v.is_consensus
                                ? 'bg-[#EAF3F5] text-[#0B4654] border border-[#B9D6DD]'
                                : 'bg-white text-[#0B4654] border border-dashed border-[#0F5C6E]'
                            }`}
                          >
                            {v.is_consensus ? 'Consensus' : 'Subclonal'}
                          </span>
                        )}
                      </div>
                    </td>
                    <td className="px-3.5 py-2.5 font-mono text-right text-[#111827] whitespace-nowrap">
                      {formatDepth(v.depth, { decimals: 0 })}
                    </td>
                    <td
                      className={`px-3.5 py-2.5 font-mono text-right font-semibold whitespace-nowrap ${
                        v.is_consensus === false ? 'text-[#4B5563]' : 'text-[#0F5C6E]'
                      }`}
                    >
                      {(v.frequency * 100).toFixed(1)}%
                    </td>
                    <td className="px-3.5 py-2.5 font-mono text-right text-[#4B5563] whitespace-nowrap">
                      {formatMetric(v.quality, { decimals: 1 })}
                    </td>
                  </tr>
                );
              })
            ) : (
              <tr>
                <td colSpan="7" className="p-8 text-center text-[#6B7280] font-medium">
                  No variants found matching the filters.
                </td>
              </tr>
            )}
          </tbody>
        </table>
      </div>

      {sortedVariants.length > 0 && (
        <div className="flex flex-col sm:flex-row sm:items-center justify-between gap-4 border-t border-[#E2E6EA] pt-3.5 text-xs text-[#6B7280] tabular-nums">
          <div className="flex flex-wrap items-center gap-4">
            <span>
              Showing <strong className="text-[#111827]">{startIndex + 1}</strong> to{' '}
              <strong className="text-[#111827]">{Math.min(startIndex + currentItemsPerPage, sortedVariants.length)}</strong> of{' '}
              <strong className="text-[#111827]">{sortedVariants.length}</strong> variants
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
                className="bg-white border border-[#CBD2D9] text-[#111827] text-[11px] rounded-[4px] px-2 py-0.5 focus:outline-none focus:border-[#0F5C6E] cursor-pointer"
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
                className="p-1.5 rounded-[4px] border border-[#CBD2D9] bg-white text-[#111827] hover:border-[#0F5C6E] disabled:opacity-40 disabled:cursor-not-allowed"
              >
                <ChevronLeft className="w-4 h-4" />
              </button>
              <span className="px-2">
                Page <strong className="text-[#111827]">{currentPage}</strong> of{' '}
                <strong className="text-[#111827]">{totalPages}</strong>
              </span>
              <button
                onClick={() => setCurrentPage((prev) => Math.min(totalPages, prev + 1))}
                disabled={currentPage === totalPages}
                className="p-1.5 rounded-[4px] border border-[#CBD2D9] bg-white text-[#111827] hover:border-[#0F5C6E] disabled:opacity-40 disabled:cursor-not-allowed"
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
