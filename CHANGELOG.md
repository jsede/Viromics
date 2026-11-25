# Changelog

All notable changes to Viromics are documented in this file.

## [2.0.0] — Milestone 1–5 Scientific & Engineering Release

### Scientific & Bioinformatics Fixes (Pipeline)
- **Decoupled Variant Calling vs. Consensus Assembly:**
  - Switched `VARIANT_CALL` to diploid likelihood calling (`bcftools call -mv --ploidy 2`), `-B` (`--no-BAQ`), `-a FORMAT/AD,FORMAT/DP,INFO/AD,INFO/ADF,INFO/ADR`, and `bcftools norm -m -any -f`.
  - Filtered reported variants using high-quality `FORMAT/DP` and per-allele `FORMAT/AD[0:1] / FORMAT/DP >= --min_freq` while decoupling consensus assembly via `prepare_consensus_vcf.py` (`AF >= --min_consensus_freq`, default `0.50`) so subclonal indels/SNPs never corrupt haploid consensus genomes.
  - Added `--use_iupac` biallelic SNP frequency gating (`[--min_iupac_freq 0.25, --max_iupac_freq 0.75)`) while keeping indels strictly major-allele.
- **Deterministic Strand-Balanced Downsampling:** Replaced coordinate-sorted `-d 8000` truncation (which caused 100% artificial strand bias on deep amplicons) with `downsample_bam_strand_balanced.py` (`--max_mpileup_depth 100000`, `--downsample_seed 42`).
- **Platform-Native Read QC (`FASTQ_QC`):** Replaced BBDuk + `fastp` with `fastplong` mean-Phred & length filtering on ONT (`-A` by default via `--ont_trim_adapters false` to preserve conserved amplicon ends) and direct raw `fastp --detect_adapter_for_pe --cut_tail --trim_poly_g` on Illumina.
- **Alignment & Primer Clipping Hardening:** Dropped unmapped, secondary, supplementary, and low-MAPQ reads in `ALIGN_READS` (`-F 0x904 -q`, `minimap2 --secondary=no`); enabled strand-aware, tolerance-windowed clipping with `samtools ampliconclip --both-ends --strand --soft-clip --tolerance --filter-len --no-excluded` plus Illumina `samtools fixmate`.
- **Deletion-Safe Depth Sync & Translation Engine (`parse_metrics.py`):**
  - Synced `samtools depth -aa -J -Q -q -G 0xF04` identically with `VARIANT_CALL`.
  - Fixed GFF3 `phase` offset subtraction on both `+` and `-` strands with partial-codon bounds checks, `-` strand adjacent frameshifts, `+1` frameshift handling, in-frame vs. frameshift indel classification, structured `consequences[]`, and staged `--cds_models` with accession aliases & slippery-site `motif_check`.
  - Implemented primer-excluded callable-target QC tiers (`PASS`/`WARN`/`FAIL`, `--qc_pass_*`/`--qc_warn_*`), two-pass primer pairing (`_LEFT_1`/`-alt1`, `_01_LEFT`, F/R), and-gated amplicon QC, reference-coordinate `low_coverage_intervals`, `compared_sites_matrices`, and `null` zero-overlap pairwise distances.
  - Added formal output schema contract [`schema/viromics_run_report.schema.json`](schema/viromics_run_report.schema.json) (`v2.0.0`) with full `run_config` and `software_versions` provenance.

### Zero-Hardcoding Dashboard (`dashboard/`)
- Added [reportSchema.js](dashboard/src/utils/reportSchema.js), [intervals.js](dashboard/src/utils/intervals.js), [phylogenetics.js](dashboard/src/utils/phylogenetics.js), and [exporters.js](dashboard/src/utils/exporters.js) eliminating hardcoded thresholds, `coverage_breadth_10x` mismatches, `150 bp` guesses, and `Date.now()` report IDs.
- Defaulted phylogenetics to **Midpoint-Rooted Neighbor-Joining** with an interactive UPGMA toggle, `compared_sites` guards, and Newick/CSV/SVG/FASTA/Multi-FASTA exports.
- Added [AmpliconChart.jsx](dashboard/src/components/AmpliconChart.jsx), indel-safe `MutationMatrix` dropout lookup (`low_coverage_intervals`), internal dropout shading (`<ReferenceArea />`), multi-track pathogen-agnostic `GeneMap`, clean print pagination, ESLint, and Vitest test coverage.

### Packaging, Schema Validation & CI
- Pinned [pipeline/environment.yml](pipeline/environment.yml) (`python=3.12`, `jsonschema=4.23.0`, `fastp=1.4.0`, `fastplong=0.7.1`, `bwa=0.7.19`, `minimap2=2.31`, `samtools=1.25`, `bcftools=1.25`, `htslib=1.25`) and removed non-functional `conda-lock.yml`.
- Added micromamba-based [Dockerfile](Dockerfile) wired to `-profile docker` and `-profile singularity` (`--container_image`), plus nf-core-compatible [nextflow_schema.json](nextflow_schema.json), offline CLI validation, and samplesheet `platform` enforcement.
- Built compact (`92 KB`) deterministic synthetic Illumina + ONT fixture (`pipeline/test_data/` + `verify_e2e.py`), behavioral unit tests, and GitHub Actions CI ([test-pipeline.yml](.github/workflows/test-pipeline.yml)).

## Roadmap — Planned for [2.1.0]
- Optional `--variant_caller clair3` (ONT haploid R10/R9) and `--variant_caller ivar` (deep Illumina subclonal calling) modules validated end-to-end under containerized Docker/Singularity CI runners.
