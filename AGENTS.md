# Viromics Developer & Agent Guidelines

## 1. Project Overview
Viromics is an integrated genomic epidemiology platform for viral amplicon whole-genome sequencing (WGS).
- **Pipeline**: Nextflow DSL2 (`pipeline/main.nf`, `pipeline/nextflow.config`, `pipeline/modules/`, `pipeline/bin/parse_metrics.py`).
- **Dashboard**: Client-side React 18 + Vite + Tailwind CSS application (`dashboard/`). Zero backend database; consumes `viromics_run_report.json`.
- **Curated Models**: `resources/annotation_models/cds_models.json` for programmed frameshifts (e.g. SARS-CoV-2 ORF1ab).
- **Interactive Wiki Documentation**: `wiki.html` (standalone interactive HTML manual in the repository root).

---

## 2. Documentation & Wiki Maintenance Policy (CRITICAL)
**Always keep `wiki.html` up to date:**
- Whenever any changes are made to:
  1. CLI arguments, flags, or defaults in `pipeline/main.nf` or `pipeline/nextflow.config`.
  2. Bioinformatics tools, versions, or process options in `pipeline/modules/`.
  3. Quality classification thresholds (PASS/WARN/FAIL), minimum coverage depth, or allele frequency filters.
  4. Curated CDS models (`resources/annotation_models/cds_models.json`) or translation exception logic.
  5. Distance matrix calculation, deletion masking, or UPGMA tree clustering logic.
  6. Output directory layout, file naming conventions, or `viromics_run_report.json` schema.
  7. Dashboard views, components, or feature tabs in `dashboard/src/`.
- **Action Required**: Immediately review and update the corresponding sections, tables, code snippets, and interactive command generator in [`wiki.html`](file:///Users/johnsebastiane/Projects/Viromics/wiki.html).

---

## 3. Scientific Invariants & Assumptions
- **Haploid Viral Genomes**: Enforce haploid variant calling (`--ploidy 1`) with optional IUPAC ambiguity codes (`--use_iupac`) for subclonal variation.
- **Fail-Closed Translation**: GFF translation exceptions without an explicitly matching curated model must return `"AA annotation unavailable"` with a distinct failure status. Never guess amino acid consequences for unsupported frameshifts.
- **Deletion Protection**: Always use `samtools depth -aa -J` so genuine biological deletions are not falsely masked with `N`.
- **SNP Distance Calculation**: Only compare coordinates covered $\ge \text{min_coverage}$ in both samples; mask structural deletions and uncalled `N`s; assign $1.0$ (100% divergence) if zero mutually covered bases exist.
- **Phylogenetic Tree Guard**: Exclude samples with $< 80\%$ coverage breadth from UPGMA tree reconstruction to prevent topological distortion. Divergence axis scale must reflect true pairwise distance ($2\times$ node height).

---

## 4. Verification Suite
Before concluding any feature or bug fix:
1. Run Python test suites:
   ```bash
   python3 pipeline/bin/test_translation_phase.py
   python3 pipeline/bin/regression_test_distance.py
   python3 pipeline/bin/test_variant_qual_filter.py
   python3 pipeline/bin/sanity_check_rsv.py
   ```
2. Build the React dashboard:
   ```bash
   npm --prefix dashboard run build
   ```
3. Verify that `wiki.html` contains up-to-date documentation.
