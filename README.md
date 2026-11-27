# Viromics: Genomic Epidemiology Suite

An integrated suite for genomic epidemiology of viral pathogens. It combines a modular Nextflow (DSL2) analytical pipeline for amplicon-based whole-genome sequencing (WGS) with a client-side React web dashboard to explore quality control (QC), genome coverage depth, variant metrics, and evolutionary relationships.

> [!TIP]
> **Interactive Documentation Wiki**: Open [`wiki.html`](wiki.html) in your browser (or run `open wiki.html` on macOS) for an interactive technical wiki with parameter builders, curated frameshift model guides, tool details, and troubleshooting FAQs!

---

## Key Features

### Nextflow DSL2 Analytical Pipeline
- **Platform-Specific Optimization**: Dedicated QC (`fastp` for paired-end Illumina with `--detect_adapter_for_pe --cut_tail --trim_poly_g`; `fastplong` mean-Phred & length filtering with `-A` by default on ONT amplicons), alignment filtering (`-F 0x904 -q`), and dynamic platform defaults.
- **Primer & Amplicon Analysis**: Strand-aware, tolerance-adjusted primer clipping (`samtools ampliconclip --both-ends --strand --soft-clip --tolerance`) with paired-end mate fixing (`samtools fixmate`), two-pass primer pairing (ARTIC `_LEFT_1`/`-alt1`, pre-direction `_01_LEFT`, F/R), and trimmed-insert callable-target QC tiers (`PASS`/`WARN`/`FAIL`).
- **Variant Calling & Translation**: Hardened `bcftools` (`-B` no-BAQ, configurable `--ont_mpileup_config`, `bcftools norm -m -any`, deterministic strand-balanced downsampling on ultra-deep amplicons), `FORMAT/AD` + `FORMAT/DP` filtering, GFF3 `phase`-accurate translation with staged `--cds_models` + slippery-site verification, in-frame vs. frameshift indel classification, and structured overlapping-gene consequences.
- **Resilient Consensus Assembly**: Strict decoupling of reported subclonal variants (`AF >= min_freq`) from consensus assembly (`AF >= min_consensus_freq`, default `0.50`), optional biallelic SNP IUPAC ambiguity coding (`--use_iupac`), and deletion-safe low-coverage masking (`samtools depth -aa -J -Q -q -G 0xF04`) with reference-coordinate `low_coverage_intervals`.
- **Epidemiological Distance & Schema Provenance**: Whole-genome SNP distance matrices (`distance_matrices`, `snp_distance_matrices`, `compared_sites_matrices`) emitting `null` for zero-overlap pairs, full `run_config` + `software_versions` provenance, parameter validation against [`nextflow_schema.json`](nextflow_schema.json), and report validation against [`schema/viromics_run_report.schema.json`](schema/viromics_run_report.schema.json).
- **Reproducible Packaging & Synthetic CI**: Pinned [`pipeline/environment.yml`](pipeline/environment.yml) (`-profile conda` / `-profile micromamba`), unified [`Dockerfile`](Dockerfile) (`-profile docker` / `-profile singularity`), and a deterministic synthetic E2E test suite (`pipeline/test_data/` + `verify_e2e.py`).

### React Web Dashboard
- **Zero-Database, Schema-Driven Design**: Client-side application that dynamically consumes `viromics_run_report.json` (`v2.0.0` schema with `v1` fallback), adapting QC tiers, thresholds, run configuration, and software provenance automatically.
- **Dynamic Multi-Strain Visualization & Amplicon QC**: Adapts to multi-FASTA reference runs, flags co-detections, shades internal coverage dropouts, and visualizes per-amplicon depth and coverage tiers.
- **Integrated Linear Genome Browser & Indel-Safe Matrix**: Multi-track staggered SVG gene map with exon connectors, interactive variant tables with structured `consequences[]`, and an indel-safe mutational sharing matrix powered by `low_coverage_intervals`.
- **Midpoint-Rooted NJ & UPGMA Phylogenetics**: Defaults to Midpoint-Rooted Neighbor-Joining with an interactive UPGMA toggle, backed by whole-genome SNP distances, pairwise `compared_sites` overlap guards, and Newick/CSV/SVG/FASTA exports.
- **Epidemiological Safe-Guard**: Excludes samples below the coverage breadth threshold (default $\ge 80\%$) or without sufficient mutually covered bases from tree reconstruction to prevent topological distortions.
- **Genomics Report Card**: A printable, multi-page paginated genomics overview summarizing deterministic run provenance, sequencing quality assurance, and Consensus vs. Subclonal mutations.

---

## Project Structure

- `pipeline/`: Modular Nextflow (DSL2) workflow script and processes.
  - `main.nf`: Core workflow script with parameter schema validation and samplesheet `platform` enforcement.
  - `nextflow.config`: Parameter sets, process resource limits, and `conda`/`micromamba`/`docker`/`singularity` profiles.
  - `environment.yml`: Pinned Conda specification (`python=3.12`, `jsonschema`, `fastp`, `fastplong`, `bwa`, `minimap2`, `samtools`, `bcftools`, `htslib`).
  - `modules/`: Decoupled Nextflow process definitions.
  - `bin/parse_metrics.py`: Metric aggregator, variant translator, and distance matrix calculator.
  - `test_data/`: Compact (`<100 KB`) deterministic Illumina + ONT synthetic fixture and `verify_e2e.py` biological assertion suite.
- `Dockerfile`: Reproducible micromamba image built from `pipeline/environment.yml` for `-profile docker` and `-profile singularity`.
- `nextflow_schema.json`: nf-core-compatible JSON Schema covering all pipeline parameters.
- `wiki.html`: Standalone interactive technical documentation wiki and manual.
- `dashboard/`: A client-side web application (React, Vite, Tailwind CSS).
- `scripts/`: Development and testing tools.
  - `generate_mock_data.py`: Synthetic multi-strain RSV-A/B dataset generator.

---

## Setup and Installation

### 1. Pre-requisites
- **Nextflow** (&ge; 23.04; tested on v26+)
- **Conda**, **Micromamba**, or **Docker** (or **Singularity**)
- **Node.js** (v18+) and **npm**

### 2. Container Execution (Docker / Singularity)
`-profile docker` and `-profile singularity` pull the prebuilt image **`ghcr.io/jsede/viromics:2.0.2`** automatically from GitHub Container Registry (GHCR) — no manual build is required. Pass `--container_image <image_or_sif>` to pin a different tag, registry image, or local Singularity `.sif` file.

For local development or custom Dockerfile edits, you can optionally build and run a local image instead:
```bash
docker build -t viromics:2.0.2 .
# Then pass: -profile docker --container_image viromics:2.0.2
```

> [!NOTE]
> If you prefer Conda or Micromamba instead of containers, pass `-profile conda` or `-profile micromamba` and Nextflow will provision the pinned environment from `pipeline/environment.yml` automatically.

### 3. Dashboard Compilation
Navigate to the `dashboard/` directory, install dependencies, and build the production bundle:
```bash
cd dashboard
npm install
npm run build
```
The optimized web application files will be written to `dashboard/dist/`.

---

## Pipeline Usage & Configuration

To run the Nextflow pipeline, execute:
```bash
nextflow run pipeline/main.nf \
  --platform [illumina|ont] \
  --reads [fastq-pattern] \
  --reference [fasta-file] \
  [optional arguments]
```

### Required Parameters
| Parameter | Description |
| :--- | :--- |
| `--platform` | Sequencing technology: `illumina` or `ont` |
| `--reads` | Path pattern to input FastQ files (e.g., `"data/*_R{1,2}.fastq.gz"` or `"data/*.fastq.gz"`) |
| `--reference` | Path to reference genome FASTA (supports multi-FASTA reference files) |

### Optional Parameters
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `--run_id` | `viromics_<sessionId>` | Deterministic run identifier for report provenance and `-resume` caching |
| `--gff` | `null` | Path to GFF3 annotation file for coding effect translation |
| `--primer_bed` | `null` | Path to primer scheme BED file for amplicon clipping and callable-target QC |
| `--amplicons_bed` | `null` | Explicit amplicon coordinates BED file (inferred from primer BED if null) |
| `--metadata` | `null` | Path to sample metadata CSV file |
| `--cds_models` | `resources/annotation_models/cds_models.json` | Curated CDS ribosomal-slippage JSON models staged explicitly to `AGGREGATE_METRICS` |
| `--min_coverage` | `10` (Illumina) / `50` (ONT) | Minimum HQ depth (`FORMAT/DP`) for variant calling |
| `--min_freq` | `0.2` (Illumina) / `0.5` (ONT) | Minimum alternative allele frequency to report variants |
| `--min_consensus_depth` | Falls back to `--min_coverage` | Minimum HQ depth to call consensus bases; lower depth is masked as `N` |
| `--min_consensus_freq` | `0.50` | Minimum allele frequency for haploid consensus & major indels |
| `--min_map_qual` | `20` | Minimum MAPQ for alignment, depth, and variant calling |
| `--min_base_qual` | `20` (Illumina) / `15` (ONT) | Minimum base quality for depth masking and variant calling |
| `--min_variant_qual` | `20` | Minimum VCF `QUAL` threshold |
| `--primer_tolerance` | `5` (Illumina) / `15` (ONT) | `samtools ampliconclip --tolerance` in bp |
| `--min_clipped_len` | `30` (Illumina) / `50` (ONT) | Minimum clipped read length retained |
| `--ont_mpileup_config` | `'ont-sup'` | `bcftools mpileup -X` profile for ONT (`'ont-sup'`, `'ont'`, or `'none'`) |
| `--ont_trim_adapters` | `false` | Enable `fastplong` de novo adapter auto-detection (disabled via `-A` by default to protect amplicon ends) |
| `--ont_min_mean_qual` | `10` | Minimum mean read Phred quality for ONT `fastplong` QC |
| `--ont_min_length` / `--ont_max_length` | `100` / `0` (unlimited) | ONT read length bounds for `fastplong` QC |
| `--ont_adapter_fasta` | `null` | Optional custom ONT adapter FASTA (skipped unless provided or `--ont_trim_adapters true`) |
| `--disable_baq` | `true` | Pass `-B` (`--no-BAQ`) to `bcftools mpileup` |
| `--max_mpileup_depth` | `100000` | Max pileup depth; triggers deterministic strand-balanced downsampling if exceeded |
| `--downsample_seed` | `42` | Fixed RNG seed for deterministic strand-balanced downsampling |
| `--use_iupac` | `false` | Enable IUPAC ambiguity codes for mixed biallelic SNPs in consensus |
| `--min_iupac_freq` / `--max_iupac_freq` | `0.25` / `0.75` | Allele frequency band for IUPAC biallelic SNP coding |
| `--qc_pass_depth` / `--qc_pass_breadth` / `--qc_pass_max_n` | `100` / `95.0` / `5.0` | Callable-target QC `PASS` thresholds (mean depth, breadth %, internal N %) |
| `--qc_warn_depth` / `--qc_warn_breadth` / `--qc_warn_max_n` | `25` / `80.0` / `10.0` | Callable-target QC `WARN` thresholds (mean depth, breadth %, internal N %) |
| `--container_image` | `'ghcr.io/jsede/viromics:2.0.2'` | Container image tag or Singularity `.sif` path for `-profile docker` / `-profile singularity` |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `'32.GB'` / `'24.h'` | Maximum per-process resource ceilings enforced across all tasks |
| `--outdir` | `./results` | Target directory to write execution results |
| `-profile` | `standard` | Environment configuration profile (`conda`, `micromamba`, `docker`, `singularity`, or `standard`) |

---

## Setting Up a Run (Step-by-Step)

### 1. Prepare Reference Files
Ensure your reference FASTA contig header names exactly match the chromosome/sequence column in your GFF3 file and the first column of your BED files:
- **Reference FASTA**: `genome.fasta`
- **Gene Annotation (GFF3)**: `annotations.gff3` containing `gene` or `CDS` features.
- **Primer Schemes (BED)**: `primers.bed` (6-column format containing primer names with `_LEFT`/`_RIGHT` or `_F`/`_R` suffixes).

### 2. Prepare Samplesheet & Metadata CSVs (Optional)
If using `--sample_sheet` instead of `--reads`, supply `sample_id,platform,fastq_1[,fastq_2]`. When the `platform` column is present (`illumina` or `ont`), it is validated against `--platform` so mixed or conflicting samplesheets fail fast. Relative FASTQ paths resolve against either the current working directory or the samplesheet directory.

To annotate samples inside the dashboard with host, collection date, or country data, create a comma-separated CSV with a mandatory `sample_id` column:
```csv
sample_id,collection_date,location,host
barcode25,2026-04-01,Melbourne,Human
barcode26,2026-04-01,Sydney,Bat
```

### 3. Run Pipeline Execution
Example execution using Oxford Nanopore multiplexed reads and a primer BED scheme under the Conda profile (substitute `-profile micromamba`, `-profile docker`, or `-profile singularity` — Docker and Singularity pull `ghcr.io/jsede/viromics:2.0.2` automatically):
```bash
nextflow run pipeline/main.nf \
  --platform ont \
  --reads "raw_data/*.fastq.gz" \
  --reference "reference/ref.fasta" \
  --gff "reference/ref.gff3" \
  --primer_bed "reference/primer_scheme.bed" \
  --metadata "metadata.csv" \
  --min_coverage 100 \
  -profile conda
```

### 4. Output Generation
The pipeline writes results for individual samples and compiles a structured run report inside the output directory:
```text
results/
  viromics_run_report.json   <-- Drag and drop this file into the React dashboard
  pipeline_info/             <-- Execution report and timeline
  [sample_id]/
    qc/                      <-- fastp quality control logs and HTML
    alignment/               <-- Sorted and primer-trimmed BAM alignment files
    variants/                <-- Filtered VCF variant sheets
    consensus/               <-- Assembled consensus FASTA genome
```

---

## Running the SARS-CoV-2 Illumina Example Dataset

The repository includes a ready-to-run SARS-CoV-2 Illumina paired-end amplicon sequencing dataset inside the `example_data/illumina_SARS2/` directory, showcasing the use of samplesheet CSV inputs and metadata annotations.

To run the pipeline on this example dataset under the Conda profile, execute:
```bash
nextflow run pipeline/main.nf \
  --platform illumina \
  --sample_sheet example_data/illumina_SARS2/samplesheet.csv \
  --reference example_data/illumina_SARS2/refs/ref.SARS-CoV-2.fasta \
  --gff example_data/illumina_SARS2/refs/ref.SARS-CoV-2.gff \
  --primer_bed example_data/illumina_SARS2/refs/primer_scheme.SARS-CoV-2-Eden.bed \
  --amplicons_bed example_data/illumina_SARS2/refs/amplicons.SARS-CoV-2-Eden.bed \
  --metadata example_data/illumina_SARS2/metadata.csv \
  --outdir example_data/illumina_SARS2/results \
  -profile conda
```

This execution:
1. Parses the samplesheet `samplesheet.csv` to map the paired-end fastq files for the negative control (`NTC`) and three clinical samples (`Sample01`, `Sample02`, `Sample03`).
2. Trims PCR primers dynamically and calls variants against the SARS-CoV-2 reference sequence.
3. Incorporates collection dates, geographic locations, and specimen host metadata from `metadata.csv`.
4. Runs the Python schema validator and outputs the consolidated run report at `example_data/illumina_SARS2/results/viromics_run_report.json`. Drag and drop this output report directly into the React dashboard to analyze the results!

---

## Running the RSV Oxford Nanopore Example Dataset

The repository also includes a scrubbed, deterministically downsampled (`~300x` mean depth per sample) Oxford Nanopore RSV amplicon dataset inside `example_data/ont_RSV/`, featuring multi-strain co-detection references (`NC_038235.1_RSV-A` and `NC_001781.1_RSV-B`):

```bash
nextflow run pipeline/main.nf \
  --platform ont \
  --sample_sheet example_data/ont_RSV/samplesheet.csv \
  --reference example_data/ont_RSV/refs/ref.RSV.fasta \
  --gff example_data/ont_RSV/refs/ref.RSV.gff \
  --primer_bed example_data/ont_RSV/refs/primer_scheme.RSV-Eden.bed \
  --amplicons_bed example_data/ont_RSV/refs/amplicons.RSV-Eden.bed \
  --metadata example_data/ont_RSV/metadata.csv \
  --outdir example_data/ont_RSV/results \
  -profile conda
```

Validate the generated report against schema `v2.0.0` and VCF decoupling invariants with:
```bash
python3 pipeline/bin/sanity_check_rsv.py example_data/ont_RSV/results
```

---

## Local Development & Testing

To test the dashboard visualizations without executing NGS runs, generate synthetic multi-strain datasets:

1. **Compile Mock Dataset**:
   ```bash
   python scripts/generate_mock_data.py
   ```
   This generates `dashboard/public/mock_run_report.json` representing a multiplexed Respiratory Syncytial Virus (RSV-A and RSV-B) run.

2. **Lint, Test, & Launch Development Server**:
   ```bash
   cd dashboard
   npm run lint
   npm run test
   npm run dev
   ```
   Open [http://localhost:5173](http://localhost:5173), upload `mock_run_report.json` (or your pipeline's `viromics_run_report.json`), and interact with the comparative heatmaps, amplicon QC charts, genome maps, and phylogenetic trees.

---

## License

Viromics is open-source software licensed under the [GNU General Public License v3.0](LICENSE) ([github.com/jsede/Viromics](https://github.com/jsede/Viromics)).
