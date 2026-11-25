# Viromics: Genomic Epidemiology Suite

An integrated suite for genomic epidemiology of viral pathogens. It combines a modular Nextflow (DSL2) analytical pipeline for amplicon-based whole-genome sequencing (WGS) with a client-side React web dashboard to explore quality control (QC), genome coverage depth, variant metrics, and evolutionary relationships.

> [!TIP]
> **Interactive Documentation Wiki**: Open [`wiki.html`](wiki.html) in your browser (or run `open wiki.html` on macOS) for an interactive technical wiki with parameter builders, curated frameshift model guides, tool details, and troubleshooting FAQs!

---

## Key Features

### Nextflow DSL2 Analytical Pipeline
- **Platform-Specific Optimization**: Dedicated QC (`fastp` for paired-end Illumina with `--detect_adapter_for_pe --cut_tail --trim_poly_g`; `fastplong` mean-Phred & length filtering for ONT), alignment filtering (`-F 0x904 -q`), and dynamic platform defaults.
- **Primer & Amplicon Analysis**: Strand-aware, tolerance-adjusted primer clipping (`samtools ampliconclip --both-ends --strand --soft-clip --tolerance`) with paired-end mate fixing (`samtools fixmate`) and amplicon coverage metrics.
- **Variant Calling**: Hardened `bcftools` (`-B` no-BAQ, `-X ont-sup` / `-X illumina`, `bcftools norm -m -any`, deterministic strand-balanced downsampling on ultra-deep amplicons), filtering on `QUAL`, high-quality `FORMAT/DP`, and per-allele `FORMAT/AD`, with strand-aware codon translation.
- **Resilient Consensus Assembly**: Strict decoupling of reported subclonal variants (`AF >= min_freq`) from consensus assembly (`AF >= min_consensus_freq`, default `0.50`), optional biallelic SNP IUPAC ambiguity coding (`--use_iupac`), and deletion-safe low-coverage masking (`samtools depth -aa -J -Q -q -G 0xF04`) synced with variant calling.
- **Epidemiological Distance Engine**: Whole-genome, pairwise genetic distance calculation. Computes SNP-only divergence across coordinates covered $\ge 10\text{x}$ in both compared samples, rendering the distance matrix frame-shift-proof and indel-resilient.

### React Web Dashboard
- **Zero-Database Design**: Client-side application that processes the pipeline's structured run report (`viromics_run_report.json`) entirely within the browser.
- **Dynamic Multi-Strain Visualization**: Automatically adapts to multi-FASTA reference runs, enabling co-detection flagging and reference strain filtering across tabs.
- **Integrated Linear Genome Browser**: Staggered gene annotation track with coordinate alignment and interactive variant hover cards.
- **Phylogenetic Tree Reconstruction**: Pre-calculated Whole-Genome SNP-only distance matrix parsed to construct a rooted UPGMA tree.
- **Epidemiological Safe-Guard**: Restricts UPGMA tree placement to samples with $\ge 80\%$ coverage breadth to prevent topological distortions from low-coverage dropouts.
- **Genomics Report Card**: A printable genomics overview summarizing sample metrics, sequencing quality assurance, and detected mutations with synonymous and non-synonymous visual indicators.

---

## Project Structure

- `pipeline/`: Modular Nextflow (DSL2) workflow script and processes.
  - `main.nf`: Core workflow script orchestrating parameters and imports.
  - `nextflow.config`: Parameter sets, resources, and conda environment profiles.
  - `modules/`: Decoupled Nextflow process definitions.
  - `bin/parse_metrics.py`: Metric aggregator, variant translator, and distance matrix calculator.
- `wiki.html`: Standalone interactive technical documentation wiki and manual.
- `dashboard/`: A client-side web application (React, Vite, Tailwind CSS).
- `scripts/`: Development and testing tools.
  - `generate_mock_data.py`: Synthetic multi-strain RSV-A/B dataset generator.

---

## Setup and Installation

### 1. Pre-requisites
- **Nextflow** (v26+)
- **Conda** or **Docker**
- **Node.js** (v18+) and **npm**

### 2. Dashboard Compilation
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
| `--gff` | `null` | Path to GFF3 annotation file for coding effect translation |
| `--primer_bed` | `null` | Path to primer scheme BED file for amplicon clipping and stats |
| `--amplicons_bed` | `null` | Explicit amplicon coordinates BED file (inferred from primer BED if null) |
| `--metadata` | `null` | Path to sample metadata CSV file |
| `--min_coverage` | `10` (Illumina) / `50` (ONT) | Minimum HQ depth (`FORMAT/DP`) for variant calling |
| `--min_freq` | `0.2` (Illumina) / `0.5` (ONT) | Minimum alternative allele frequency to report variants |
| `--min_consensus_depth` | Falls back to `--min_coverage` | Minimum HQ depth to call consensus bases; lower depth is masked as `N` |
| `--min_consensus_freq` | `0.50` | Minimum allele frequency for haploid consensus & major indels |
| `--min_map_qual` | `20` | Minimum MAPQ for alignment, depth, and variant calling |
| `--min_base_qual` | `20` (Illumina) / `15` (ONT) | Minimum base quality for depth masking and variant calling |
| `--min_variant_qual` | `20` | Minimum VCF `QUAL` threshold |
| `--primer_tolerance` | `5` (Illumina) / `15` (ONT) | `samtools ampliconclip --tolerance` in bp |
| `--min_clipped_len` | `30` (Illumina) / `50` (ONT) | Minimum clipped read length retained |
| `--ont_min_mean_qual` | `10` | Minimum mean read Phred quality for ONT `fastplong` QC |
| `--ont_min_length` / `--ont_max_length` | `100` / `0` (unlimited) | ONT read length bounds for `fastplong` QC |
| `--ont_adapter_fasta` | `null` | Optional custom ONT adapter FASTA (auto-detected if null) |
| `--disable_baq` | `true` | Pass `-B` (`--no-BAQ`) to `bcftools mpileup` |
| `--max_mpileup_depth` | `100000` | Max pileup depth; triggers deterministic strand-balanced downsampling if exceeded |
| `--downsample_seed` | `42` | Fixed RNG seed for deterministic strand-balanced downsampling |
| `--use_iupac` | `false` | Enable IUPAC ambiguity codes for mixed biallelic SNPs in consensus |
| `--min_iupac_freq` / `--max_iupac_freq` | `0.25` / `0.75` | Allele frequency band for IUPAC biallelic SNP coding |
| `--outdir` | `./results` | Target directory to write execution results |
| `-profile` | `standard` | Environment configuration profile (e.g., `-profile conda`) |

---

## Setting Up a Run (Step-by-Step)

### 1. Prepare Reference Files
Ensure your reference FASTA contig header names exactly match the chromosome/sequence column in your GFF3 file and the first column of your BED files:
- **Reference FASTA**: `genome.fasta`
- **Gene Annotation (GFF3)**: `annotations.gff3` containing `gene` or `CDS` features.
- **Primer Schemes (BED)**: `primers.bed` (6-column format containing primer names with `_LEFT`/`_RIGHT` or `_F`/`_R` suffixes).

### 2. Prepare Metadata CSV (Optional)
To annotate samples inside the dashboard with host, collection date, or country data, create a comma-separated CSV with a mandatory `sample_id` column:
```csv
sample_id,collection_date,location,host
barcode25,2026-04-01,Melbourne,Human
barcode26,2026-04-01,Sydney,Bat
```

### 3. Run Pipeline Execution
Example execution using Oxford Nanopore multiplexed reads and a primer BED scheme under the Conda profile:
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

## Local Development & Testing

To test the dashboard visualizations without executing NGS runs, generate synthetic multi-strain datasets:

1. **Compile Mock Dataset**:
   ```bash
   python scripts/generate_mock_data.py
   ```
   This generates `dashboard/public/mock_run_report.json` representing a multiplexed Respiratory Syncytial Virus (RSV-A and RSV-B) run.

2. **Launch Development Server**:
   ```bash
   cd dashboard
   npm run dev
   ```
   Open [http://localhost:5173](http://localhost:5173), upload `mock_run_report.json`, and interact with the comparative heatmaps, alignments, and trees.
