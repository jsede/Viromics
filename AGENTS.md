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
  5. Distance matrix calculation, deletion masking, or NJ/UPGMA tree clustering logic.
  6. Output directory layout, file naming conventions, or `viromics_run_report.json` schema.
  7. Dashboard views, components, or feature tabs in `dashboard/src/`.
- **Action Required**: Immediately review and update the corresponding sections, tables, code snippets, and interactive command generator in [`wiki.html`](file:///Users/johnsebastiane/Projects/Viromics/wiki.html).

---

## 3. Scientific Invariants & Assumptions
- **Haploid Viral Genomes**: Perform diploid likelihood variant calling (`bcftools call -mv --ploidy 2` + `bcftools norm -m -any`) to retain true subclonal alleles (`AF >= --min_freq`), while strictly decoupling consensus assembly (`AF >= --min_consensus_freq`, default `0.50`) and optional biallelic SNP IUPAC ambiguity calling (`--use_iupac`).
- **Fail-Closed Translation**: GFF3 `phase` offsets are subtracted from the coding 5' end on both `+` and `-` strands with strict partial-codon bounds checks. GFF translation exceptions without an explicitly matching curated model (including reference alias & motif verification) must return `"AA annotation unavailable"` with a distinct failure status.
- **Deletion Protection & Filter Sync**: Always run `samtools depth -aa -J -Q <min_map_qual> -q <min_base_qual> -G 0xF04` synced identically with `VARIANT_CALL` so genuine biological deletions are not falsely masked with `N` and low-MAPQ/BaseQ reads never inflate consensus breadth.
- **Callable-Target QC Tiers**: When `--primer_bed` or `--amplicons_bed` is supplied, evaluate QC tiers (`PASS`/`WARN`/`FAIL`, configurable via `--qc_pass_*` / `--qc_warn_*`) on the union of primer-trimmed insert intervals (`[left_end, right_start]`) and internal `N%`, while reporting both whole-contig and callable-target metrics.
- **SNP Distance Calculation**: Only compare coordinates covered $\ge \text{min_consensus_depth}$ in both samples; mask structural deletions and uncalled `N`s; emit `null` (not `1.0`) alongside `compared_sites_matrices = 0` when zero mutually covered bases exist.
- **Phylogenetic Tree Guard**: Exclude samples below the coverage breadth threshold ($< 80\%$ default) or with `null` / insufficient mutually covered sites from tree reconstruction to prevent topological distortion. Default to Midpoint-Rooted Neighbor-Joining with an optional UPGMA toggle; divergence axis scale must reflect true pairwise distance.

---

## 4. Verification Suite
Before concluding any feature or bug fix:
1. Run Python test suites:
   ```bash
   python3 pipeline/bin/test_prepare_consensus_vcf.py
   python3 pipeline/bin/test_translation_phase.py
   python3 pipeline/bin/regression_test_distance.py
   python3 pipeline/bin/test_samtools_aa.py
   python3 pipeline/bin/test_variant_qual_filter.py
   python3 pipeline/bin/test_parse_metrics_m2.py
   ```
2. Lint, test, and build the React dashboard:
   ```bash
   npm --prefix dashboard run lint
   npm --prefix dashboard run test
   npm --prefix dashboard run build
   ```
3. Verify that `wiki.html` contains up-to-date documentation, `nextflow_schema.json` matches `pipeline/nextflow.config`, `schema/viromics_run_report.schema.json` validates report outputs, and synthetic E2E verification (`pipeline/test_data/verify_e2e.py`) passes when workflow or metrics logic changes.
