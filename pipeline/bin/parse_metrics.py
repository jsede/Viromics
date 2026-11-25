#!/usr/bin/env python3
"""
Viromics metric aggregator, variant consequence annotator, callable-target QC
classifier, pairwise SNP distance engine, and schema validator.
"""

import argparse
import csv
import json
import os
import re
import subprocess
import sys
import urllib.parse
from collections import defaultdict

SCHEMA_VERSION = "2.0.0"

CURATED_CDS_MODELS_CACHE = None
CURATED_CDS_MODELS_PATH = None


def load_curated_models_cached(cds_models_path=None):
    """Load curated CDS frameshift models.

    Explicitly staged paths fail loudly if missing or invalid JSON.
    When invoked without an explicit path (e.g. unit tests), falls back to
    the bundled repository models and respects CURATED_CDS_MODELS_CACHE.
    """
    global CURATED_CDS_MODELS_CACHE

    if cds_models_path is None and CURATED_CDS_MODELS_PATH is None:
        if CURATED_CDS_MODELS_CACHE is not None:
            return CURATED_CDS_MODELS_CACHE

    explicit_path = cds_models_path if cds_models_path is not None else CURATED_CDS_MODELS_PATH
    if explicit_path:
        if not os.path.exists(explicit_path):
            raise FileNotFoundError(f"Curated CDS models file not found: {explicit_path}")
        try:
            with open(explicit_path, "r") as f:
                models = json.load(f)
        except Exception as exc:
            raise ValueError(f"Failed to parse curated CDS models JSON '{explicit_path}': {exc}") from exc
        if not isinstance(models, list):
            raise ValueError(f"Curated CDS models JSON '{explicit_path}' must contain a JSON array.")
        if cds_models_path is None or cds_models_path == CURATED_CDS_MODELS_PATH:
            CURATED_CDS_MODELS_CACHE = models
        return models

    script_dir = os.path.dirname(os.path.abspath(__file__))
    paths = [
        os.path.join(script_dir, "resources", "annotation_models", "cds_models.json"),
        os.path.join(script_dir, "..", "..", "resources", "annotation_models", "cds_models.json"),
        "resources/annotation_models/cds_models.json",
    ]
    for p in paths:
        if os.path.exists(p):
            try:
                with open(p, "r") as f:
                    models = json.load(f)
            except Exception as exc:
                raise ValueError(f"Failed to parse bundled curated CDS models '{p}': {exc}") from exc
            if not isinstance(models, list):
                raise ValueError(f"Bundled curated CDS models '{p}' must contain a JSON array.")
            CURATED_CDS_MODELS_CACHE = models
            return CURATED_CDS_MODELS_CACHE

    CURATED_CDS_MODELS_CACHE = []
    return CURATED_CDS_MODELS_CACHE


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate NGS run metrics for Viromics Dashboard")
    parser.add_argument("--run-id", required=True, help="Unique identifier for the run")
    parser.add_argument("--fastp", required=True, help="Comma-separated list of fastp/fastplong JSON files")
    parser.add_argument("--bams", required=True, help="Comma-separated list of BAM files")
    parser.add_argument("--vcfs", required=True, help="Comma-separated list of VCF files")
    parser.add_argument("--consensuses", required=True, help="Comma-separated list of consensus FASTA files")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--gff", help="Optional GFF3 annotation file")
    parser.add_argument("--primers", help="Optional path to primer scheme BED file")
    parser.add_argument("--amplicons-bed", help="Optional path to explicit amplicon coordinates BED file")
    parser.add_argument("--metadata", help="Optional path to sample metadata CSV file")
    parser.add_argument("--cds-models", help="Optional explicit path to curated CDS frameshift JSON models")
    parser.add_argument("--versions", help="Comma-separated list of tool versions YAML files")
    parser.add_argument("--run-config", help="Optional path to resolved Nextflow run_config JSON file")
    parser.add_argument("--min-coverage", type=int, default=10, help="Minimum HQ depth for variant calling & coverage breadth")
    parser.add_argument("--min-freq", type=float, default=0.2, help="Minimum alternative allele frequency for variants")
    parser.add_argument("--min-consensus-depth", type=int, default=None, help="Minimum HQ depth for consensus masking (defaults to --min-coverage)")
    parser.add_argument("--min-consensus-freq", type=float, default=0.50, help="Minimum AF for major/consensus alleles")
    parser.add_argument("--min-map-qual", type=int, default=20, help="Minimum mapping quality (MAPQ) synced with consensus & variant calling")
    parser.add_argument("--min-base-qual", type=int, default=20, help="Minimum base quality (BaseQ) synced with consensus & variant calling")
    parser.add_argument("--qc-pass-depth", type=int, default=100, help="Minimum callable mean depth for QC PASS")
    parser.add_argument("--qc-pass-breadth", type=float, default=95.0, help="Minimum callable breadth (percent) at min_coverage for QC PASS")
    parser.add_argument("--qc-pass-max-n", type=float, default=5.0, help="Maximum internal N rate (percent) for QC PASS")
    parser.add_argument("--qc-warn-depth", type=int, default=25, help="Minimum callable mean depth for QC WARN")
    parser.add_argument("--qc-warn-breadth", type=float, default=80.0, help="Minimum callable breadth (percent) at min_coverage for QC WARN")
    parser.add_argument("--qc-warn-max-n", type=float, default=10.0, help="Maximum internal N rate (percent) for QC WARN")
    parser.add_argument("--output", default="viromics_run_report.json", help="Output JSON path")
    return parser.parse_args()


def parse_reference(ref_path):
    """Parse all reference sequences from a multi-FASTA file."""
    references = []
    current_name = None
    current_seq = []
    with open(ref_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_name:
                    seq = "".join(current_seq)
                    references.append({
                        "name": current_name,
                        "length": len(seq),
                        "sequence": seq,
                    })
                current_name = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_name:
            seq = "".join(current_seq)
            references.append({
                "name": current_name,
                "length": len(seq),
                "sequence": seq,
            })
    return references


def parse_gff(gff_path):
    """Parse genes and CDS from GFF3 file, grouped by reference sequence ID."""
    features = defaultdict(list)
    if not gff_path or not os.path.exists(gff_path):
        return dict(features)

    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue

            ref_id = parts[0]
            feature_type = parts[2]
            if feature_type in ["gene", "CDS", "exon"]:
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]

                phase = 0
                if feature_type == "CDS" and len(parts) >= 8:
                    phase_col = parts[7].strip()
                    if phase_col in [".", "", "null"]:
                        phase = 0
                    elif phase_col in ["0", "1", "2"]:
                        phase = int(phase_col)
                    else:
                        phase = -1  # Invalid phase

                attrs = {}
                for attr in parts[8].split(";"):
                    if not attr.strip():
                        continue
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        attrs[k.strip()] = urllib.parse.unquote(v.strip())

                name = attrs.get("gene", attrs.get("Name", attrs.get("ID", "Unknown")))
                product = attrs.get("product", "")

                features[ref_id].append({
                    "ref_id": ref_id,
                    "name": name,
                    "type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "phase": phase,
                    "product": product,
                    "id": attrs.get("ID", ""),
                    "parent": attrs.get("Parent", ""),
                    "gene": attrs.get("gene", ""),
                    "feature_name": attrs.get("Name", ""),
                    "locus_tag": attrs.get("locus_tag", ""),
                    "protein_id": attrs.get("protein_id", ""),
                    "transcript_id": attrs.get("transcript_id", ""),
                    "exception": attrs.get("exception", ""),
                    "note": attrs.get("Note", ""),
                    "transl_except": attrs.get("transl_except", ""),
                    "attributes": attrs,
                })

    # Hierarchical filtering per reference: CDS > gene > exon to prevent duplication clutter
    sorted_genes = {}
    for ref_id, f_list in features.items():
        cds_features = [f for f in f_list if f["type"] == "CDS"]
        if cds_features:
            selected = cds_features
        else:
            gene_features = [f for f in f_list if f["type"] == "gene"]
            selected = gene_features if gene_features else f_list

        cds_only = [f for f in selected if f["type"] == "CDS"]
        non_cds_only = [f for f in selected if f["type"] != "CDS"]

        cds_groups = defaultdict(list)
        ungrouped_count = 0
        for f in cds_only:
            attrs = f["attributes"]
            grouping_key = None
            for key in ["protein_id", "transcript_id", "ID", "Parent", "locus_tag"]:
                val = attrs.get(key, "").strip()
                if val and val.lower() != "unknown":
                    grouping_key = val
                    break

            if not grouping_key:
                ungrouped_count += 1
                grouping_key = f"__ungrouped_cds_{ungrouped_count}__"

            f["grouping_key"] = grouping_key
            cds_groups[(f["strand"], grouping_key)].append(f)

        groups_list = []
        for (strand, _gkey), group_feats in cds_groups.items():
            # Sort CDS parts in biological coding order:
            # + strand: ascending start, then ascending end
            # - strand: descending start, then descending end
            if strand == "+":
                group_feats.sort(key=lambda x: (x["start"], x["end"]))
            else:
                group_feats.sort(key=lambda x: (-x["start"], -x["end"]))

            distinct_intervals = {(f["start"], f["end"]) for f in group_feats}
            is_multi = len(group_feats) >= 2 and len(distinct_intervals) >= 2

            has_fs = False
            has_te = False
            for f in group_feats:
                exc_val = " ".join(f.get("exception", "").lower().replace("_", " ").split())
                note_val = " ".join(f.get("note", "").lower().replace("_", " ").split())
                if (
                    "programmed ribosomal frameshift" in exc_val
                    or "ribosomal slippage" in exc_val
                    or "frameshift" in exc_val
                    or "programmed ribosomal frameshift" in note_val
                    or "ribosomal slippage" in note_val
                    or "frameshift" in note_val
                ):
                    has_fs = True
                if f.get("transl_except", "").strip():
                    has_te = True

            if has_fs:
                exception_status = "unsupported_programmed_frameshift"
            elif has_te:
                exception_status = "unsupported_translation_exception"
            else:
                exception_status = None

            requires_curated = exception_status == "unsupported_programmed_frameshift"

            for f in group_feats:
                f["multi_part"] = is_multi
                f["translation_exception_status"] = exception_status
                f["requires_curated_model"] = requires_curated

            groups_list.append(group_feats)

        for f in non_cds_only:
            f["multi_part"] = False
            f["translation_exception_status"] = None
            f["requires_curated_model"] = False
            f["grouping_key"] = None
            groups_list.append([f])

        groups_list.sort(key=lambda g: min(f["start"] for f in g))
        sorted_genes[ref_id] = [f for g in groups_list for f in g]

    return sorted_genes


def parse_fastp(fastp_path):
    """Parse quality stats from fastp (Illumina) or fastplong (ONT) JSON report."""
    if not os.path.exists(fastp_path):
        raise FileNotFoundError(f"QC JSON report not found: {fastp_path}")

    try:
        with open(fastp_path, "r") as f:
            data = json.load(f)
    except Exception as exc:
        raise ValueError(f"Failed to parse QC JSON '{fastp_path}': {exc}") from exc

    if not isinstance(data, dict) or "summary" not in data:
        raise ValueError(f"QC JSON '{fastp_path}' is missing required 'summary' block.")

    summary = data.get("summary", {}) or {}
    before = summary.get("before_filtering", {}) or {}
    after = summary.get("after_filtering", {}) or {}

    def _opt_float(val, ndigits=6):
        if val is None:
            return None
        return round(float(val), ndigits)

    # fastp emits duplication at top-level data["duplication"]["rate"];
    # fastplong omits duplication completely -> preserve as None (null in JSON).
    dup_rate = None
    top_dup = data.get("duplication")
    if isinstance(top_dup, dict) and top_dup.get("rate") is not None:
        dup_rate = round(float(top_dup["rate"]), 6)
    else:
        sum_dup = summary.get("duplication")
        if isinstance(sum_dup, dict) and sum_dup.get("rate") is not None:
            dup_rate = round(float(sum_dup["rate"]), 6)

    q30_raw = _opt_float(before.get("q30_rate"))
    q30_trimmed = _opt_float(after.get("q30_rate"))

    mean_len_raw = before.get("read_mean_length", before.get("read1_mean_length"))
    mean_len_trimmed = after.get("read_mean_length", after.get("read1_mean_length"))

    return {
        "raw_reads": int(before.get("total_reads", 0)),
        "raw_bases": int(before.get("total_bases", 0)),
        "trimmed_reads": int(after.get("total_reads", 0)),
        "trimmed_bases": int(after.get("total_bases", 0)),
        "q30_rate_raw": q30_raw,
        "q30_rate_trimmed": q30_trimmed,
        "duplication_rate": dup_rate,
        "mean_read_length_raw": _opt_float(mean_len_raw, 2),
        "mean_read_length_trimmed": _opt_float(mean_len_trimmed, 2),
    }


def compute_low_coverage_intervals(ref_depths, threshold):
    """Return 1-based inclusive [start, end] intervals where depth < threshold."""
    intervals = []
    in_run = False
    run_start = None
    for idx, depth in enumerate(ref_depths):
        if depth < threshold:
            if not in_run:
                run_start = idx + 1
                in_run = True
        else:
            if in_run:
                intervals.append([run_start, idx])
                in_run = False
    if in_run:
        intervals.append([run_start, len(ref_depths)])
    return intervals


def calculate_coverage(
    bam_path,
    ref_lengths,
    num_bins=500,
    min_cov=10,
    min_consensus_depth=None,
    min_map_qual=20,
    min_base_qual=20,
):
    """Calculate depth of coverage synced with CONSENSUS_GEN/VARIANT_CALL filters."""
    results = {}
    mask_depth = min_consensus_depth if min_consensus_depth is not None else min_cov

    depths = {ref_name: [0] * ref_len for ref_name, ref_len in ref_lengths.items()}
    mapped_reads = {ref_name: 0 for ref_name in ref_lengths}

    # Count primary alignments passing MAPQ filter per contig (exclude unmapped, secondary, supplementary: 0x904 = 2308)
    for ref_name in ref_lengths:
        cmd = ["samtools", "view", "-c", "-F", "2308", "-q", str(min_map_qual), bam_path, ref_name]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        mapped_reads[ref_name] = int(result.stdout.strip() or "0")

    # Depth synced with CONSENSUS_GEN: -aa (all positions), -J (count deletion-spanning reads), -Q MAPQ, -q BaseQ, -G 0xF04
    depth_cmd = [
        "samtools",
        "depth",
        "-aa",
        "-J",
        "-Q",
        str(min_map_qual),
        "-q",
        str(min_base_qual),
        "-G",
        "0xF04",
        bam_path,
    ]
    depth_result = subprocess.run(depth_cmd, capture_output=True, text=True, check=True)
    for line in depth_result.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            ref_name = parts[0]
            pos = int(parts[1]) - 1
            depth = int(parts[2])
            if ref_name in depths:
                ref_len = ref_lengths[ref_name]
                if 0 <= pos < ref_len:
                    depths[ref_name][pos] = depth

    for ref_name, ref_len in ref_lengths.items():
        ref_depths = depths[ref_name]
        mean_depth = sum(ref_depths) / ref_len if ref_len > 0 else 0.0
        cov_10x = sum(1 for d in ref_depths if d >= 10) / ref_len if ref_len > 0 else 0.0
        cov_30x = sum(1 for d in ref_depths if d >= 30) / ref_len if ref_len > 0 else 0.0
        cov_min_cov = sum(1 for d in ref_depths if d >= min_cov) / ref_len if ref_len > 0 else 0.0

        bin_size = max(1, ref_len // num_bins)
        downsampled_depths = []
        for i in range(0, ref_len, bin_size):
            chunk = ref_depths[i : i + bin_size]
            avg_depth = round(sum(chunk) / len(chunk)) if chunk else 0
            downsampled_depths.append({
                "pos": i + 1,
                "depth": avg_depth,
            })

        results[ref_name] = {
            "mean_depth": round(mean_depth, 2),
            "coverage_breadth_10x": round(cov_10x * 100, 2),
            "coverage_breadth_30x": round(cov_30x * 100, 2),
            "coverage_breadth_min_cov": round(cov_min_cov * 100, 2),
            "mapped_reads": mapped_reads[ref_name],
            "coverage_profile": downsampled_depths,
            "low_coverage_intervals": compute_low_coverage_intervals(ref_depths, mask_depth),
        }
    return results, depths


GENETIC_CODE = {
    "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
    "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
    "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
}

COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
    "N": "N", "n": "n",
}


def _classify_snp_effect(wt_aa, mut_aa):
    if wt_aa == mut_aa:
        return "synonymous"
    if mut_aa == "*":
        return "nonsense"
    if wt_aa == "*":
        return "stop_lost"
    return "missense"


def find_matching_model(ref_id, gene, models):
    """Find a curated CDS frameshift model matching reference, strand, feature identifiers, and exception status."""
    if not models:
        return None
    ref_str = ref_id or ""
    for model in models:
        ref_exact = model.get("reference")
        ref_aliases = model.get("reference_aliases", []) or []
        ref_patterns = model.get("reference_patterns", []) or []

        if ref_exact or ref_aliases or ref_patterns:
            matched_ref = False
            if ref_exact and ref_str == ref_exact:
                matched_ref = True
            elif ref_str in ref_aliases:
                matched_ref = True
            else:
                for pat in ref_patterns:
                    if re.search(pat, ref_str):
                        matched_ref = True
                        break
            if not matched_ref:
                continue

        model_strand = model.get("strand")
        if model_strand and model_strand != gene.get("strand"):
            continue

        feature_names = model.get("feature_names", [])
        products = model.get("products", [])
        prot_id = gene.get("protein_id", "")
        candidate_names = {
            gene.get("name"),
            gene.get("id"),
            gene.get("gene"),
            gene.get("feature_name"),
            prot_id,
            f"cds-{prot_id}" if prot_id else None,
        }
        name_match = any(n in feature_names for n in candidate_names if n) or (gene.get("product") in products)
        if not name_match:
            continue

        model_status = model.get("translation_exception_status")
        if model_status and model_status != gene.get("translation_exception_status"):
            continue

        return model
    return None


def run_curated_aa_annotation(pos, ref, alt, model, ref_seq):
    """Translate a SNP across curated CDS intervals supporting programmed -1/+1 frameshifts on + and - strands."""
    if not ref_seq:
        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_reference_mismatch"}

    motif = model.get("motif_check")
    if isinstance(motif, dict):
        m_start = motif.get("start")
        m_end = motif.get("end")
        m_seq = motif.get("sequence", "")
        if (
            m_start is None
            or m_end is None
            or m_start < 1
            or m_end > len(ref_seq)
            or m_start > m_end
            or ref_seq[m_start - 1 : m_end].upper() != m_seq.upper()
        ):
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_reference_mismatch"}

    strand = model.get("strand", "+")
    cds_bases = []
    cds_positions = []
    prev_start = None
    prev_end = None

    for intv in model.get("intervals", []):
        start = intv.get("start")
        end = intv.get("end")
        phase = intv.get("phase", 0)
        frame_offset = intv.get("frame_offset", 0)

        if start is None or end is None or start < 1 or end > len(ref_seq) or start > end:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_reference_mismatch"}

        if frame_offset not in [-1, 0, 1]:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}

        if strand == "+":
            intv_pos = list(range(start, end + 1))
            intv_seq = list(ref_seq[start - 1 : end])
            if prev_end is not None and frame_offset == -1:
                if start == prev_end + 1:
                    # Adjacent intervals with programmed -1 frameshift: re-read prev_end base
                    intv_pos = [prev_end] + intv_pos
                    intv_seq = [ref_seq[prev_end - 1]] + intv_seq
                elif start > prev_end + 1:
                    return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}
            elif prev_end is not None and frame_offset == 1:
                if start == prev_end + 1:
                    # Adjacent intervals with programmed +1 frameshift: skip 1 base at interval start
                    if not intv_seq:
                        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}
                    intv_pos = intv_pos[1:]
                    intv_seq = intv_seq[1:]
                elif start == prev_end + 2:
                    # Explicit 1 bp genomic gap between intervals encodes +1 shift cleanly
                    pass
                else:
                    return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}
        else:
            intv_pos = list(range(end, start - 1, -1))
            raw_seq = ref_seq[start - 1 : end]
            intv_seq = [COMPLEMENT.get(b, b) for b in reversed(raw_seq)]
            if prev_start is not None and frame_offset == -1:
                if end == prev_start - 1:
                    # Adjacent - strand intervals in descending coding order: re-read prev_start base
                    intv_pos = [prev_start] + intv_pos
                    intv_seq = [COMPLEMENT.get(ref_seq[prev_start - 1], "N")] + intv_seq
                elif end < prev_start - 1:
                    return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}
            elif prev_start is not None and frame_offset == 1:
                if end == prev_start - 1:
                    if not intv_seq:
                        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}
                    intv_pos = intv_pos[1:]
                    intv_seq = intv_seq[1:]
                elif end == prev_start - 2:
                    pass
                else:
                    return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}

        prev_start = start
        prev_end = end

        if phase < 0 or phase > len(intv_seq):
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}

        intv_seq = intv_seq[phase:]
        intv_pos = intv_pos[phase:]

        cds_bases.extend(intv_seq)
        cds_positions.extend(intv_pos)

    cds_sequence = "".join(cds_bases)

    matching_indices = [i for i, p in enumerate(cds_positions) if p == pos]
    if not matching_indices:
        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_position_out_of_bounds"}

    aa_changes = []
    effects = []
    mut_base = alt if strand == "+" else COMPLEMENT.get(alt, alt)

    for cds_idx in matching_indices:
        codon_start_idx = (cds_idx // 3) * 3
        codon_offset = cds_idx % 3

        if codon_start_idx < 0 or codon_start_idx + 3 > len(cds_sequence):
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_codon_incomplete"}

        wt_codon = cds_sequence[codon_start_idx : codon_start_idx + 3]
        if len(wt_codon) != 3:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_codon_incomplete"}

        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = mut_base
        mut_codon = "".join(mut_codon_list)

        wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
        mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
        aa_pos = (cds_idx // 3) + 1

        if wt_aa == "?" or mut_aa == "?":
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}

        effects.append(_classify_snp_effect(wt_aa, mut_aa))
        if wt_aa == mut_aa:
            aa_changes.append(f"{wt_aa}{aa_pos}{mut_aa} (syn)")
        else:
            aa_changes.append(f"{wt_aa}{aa_pos}{mut_aa}")

    if not aa_changes:
        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "curated_model_boundary_unsupported"}

    return {"aa_change": "/".join(aa_changes), "effect": "/".join(dict.fromkeys(effects)), "status": None}


def get_aa_change(pos, ref, alt, gene, ref_seq, curated_models=None):
    """Calculate codon-level amino acid change for SNPs and classify coding indels."""
    if len(ref) != len(alt):
        delta = len(alt) - len(ref)
        abs_delta = abs(delta)
        is_inframe = (abs_delta % 3) == 0
        if delta < 0:
            if is_inframe:
                aa_count = abs_delta // 3
                return {"aa_change": f"In-frame del ({aa_count} aa)", "effect": "inframe_deletion", "status": None}
            return {"aa_change": f"Frameshift del (-{abs_delta} bp)", "effect": "frameshift_deletion", "status": None}
        if is_inframe:
            aa_count = abs_delta // 3
            return {"aa_change": f"In-frame ins ({aa_count} aa)", "effect": "inframe_insertion", "status": None}
        return {"aa_change": f"Frameshift ins (+{abs_delta} bp)", "effect": "frameshift_insertion", "status": None}

    if len(ref) > 1:
        return {"aa_change": f"MNP ({len(ref)} bp)", "effect": "mnp", "status": None}

    if not ref_seq:
        return {"aa_change": "", "effect": "unknown", "status": None}

    if gene.get("translation_exception_status"):
        models = curated_models if curated_models is not None else load_curated_models_cached()
        matching_model = find_matching_model(gene.get("ref_id"), gene, models)
        if not matching_model:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": gene["translation_exception_status"]}
        return run_curated_aa_annotation(pos, ref, alt, matching_model, ref_seq)

    phase = gene.get("phase", 0)
    if phase == -1:
        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "invalid_phase"}
    if gene.get("multi_part", False):
        return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "unsupported_multipart_cds"}

    g_start = gene["start"]
    g_end = gene["end"]
    strand = gene.get("strand", "+")

    if not (g_start <= pos <= g_end):
        return {"aa_change": "", "effect": "unknown", "status": None}

    if strand == "+":
        # GFF3 spec: phase (0,1,2) is the number of bases to skip at the 5' start of the CDS feature
        rel_offset = (pos - g_start) - phase
        if rel_offset < 0:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "incomplete_partial_codon"}
        codon_index = rel_offset // 3
        codon_offset = rel_offset % 3
        codon_start_1b = g_start + phase + codon_index * 3
        codon_end_1b = codon_start_1b + 2

        if codon_start_1b < g_start or codon_end_1b > g_end or codon_start_1b < 1 or codon_end_1b > len(ref_seq):
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "incomplete_partial_codon"}

        wt_codon = ref_seq[codon_start_1b - 1 : codon_end_1b]
        if len(wt_codon) != 3:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "incomplete_partial_codon"}

        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = alt
        mut_codon = "".join(mut_codon_list)
    else:
        # Negative-strand GFF3: CDS starts at g_end (5' coding) and proceeds toward g_start
        rel_offset = (g_end - pos) - phase
        if rel_offset < 0:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "incomplete_partial_codon"}
        codon_index = rel_offset // 3
        codon_offset = rel_offset % 3
        codon_hi_1b = g_end - phase - codon_index * 3
        codon_lo_1b = codon_hi_1b - 2

        if codon_lo_1b < g_start or codon_hi_1b > g_end or codon_lo_1b < 1 or codon_hi_1b > len(ref_seq):
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "incomplete_partial_codon"}

        genomic_codon = ref_seq[codon_lo_1b - 1 : codon_hi_1b]
        if len(genomic_codon) != 3:
            return {"aa_change": "AA annotation unavailable", "effect": "unknown", "status": "incomplete_partial_codon"}

        wt_codon = "".join(COMPLEMENT.get(b, b) for b in reversed(genomic_codon))
        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = COMPLEMENT.get(alt, alt)
        mut_codon = "".join(mut_codon_list)

    wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
    mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
    aa_pos = codon_index + 1

    if wt_aa == "?" or mut_aa == "?":
        return {"aa_change": "", "effect": "unknown", "status": None}

    effect = _classify_snp_effect(wt_aa, mut_aa)
    if wt_aa == mut_aa:
        return {"aa_change": f"{wt_aa}{aa_pos}{mut_aa} (syn)", "effect": effect, "status": None}
    return {"aa_change": f"{wt_aa}{aa_pos}{mut_aa}", "effect": effect, "status": None}


def _extract_vcf_dp_and_af(parts, info):
    """Extract (dp, af) from FORMAT/AD+DP first, then INFO/DP4, then INFO/DP+VAF/AF. Never fabricates 100 or 1.0."""
    if len(parts) >= 10 and parts[8]:
        fmt_keys = parts[8].split(":")
        smp_vals = parts[9].split(":")
        fmt_map = dict(zip(fmt_keys, smp_vals))
        ad_str = fmt_map.get("AD")
        dp_str = fmt_map.get("DP")
        if ad_str and ad_str != ".":
            try:
                ad_vals = [int(x) for x in ad_str.split(",")]
                if len(ad_vals) >= 2:
                    alt_ad = ad_vals[1]
                    if dp_str and dp_str != ".":
                        dp = int(dp_str)
                    else:
                        dp = sum(ad_vals)
                    af = (alt_ad / dp) if dp > 0 else 0.0
                    return dp, af
            except ValueError:
                pass

    if "DP4" in info and isinstance(info["DP4"], str):
        try:
            dp4 = [int(x) for x in info["DP4"].split(",")]
            if len(dp4) == 4:
                ref_cnt = dp4[0] + dp4[1]
                alt_cnt = dp4[2] + dp4[3]
                dp = ref_cnt + alt_cnt
                af = (alt_cnt / dp) if dp > 0 else 0.0
                return dp, af
        except ValueError:
            pass

    dp = 0
    if "DP" in info and isinstance(info["DP"], str) and info["DP"] != ".":
        try:
            dp = int(info["DP"])
        except ValueError:
            dp = 0

    af = 0.0
    for af_key in ("VAF", "AF"):
        if af_key in info and isinstance(info[af_key], str) and info[af_key] != ".":
            try:
                af = float(info[af_key].split(",")[0])
                break
            except ValueError:
                pass

    return dp, af


def _get_feature_display_name(gene):
    """Disambiguate shared gene=ORF1ab vs product=ORF1a polyprotein labels cleanly."""
    name = gene.get("name", "Unknown")
    product = (gene.get("product") or "").lower()
    if name == "ORF1ab" and "orf1a polyprotein" in product:
        return "ORF1a"
    return name


def parse_vcf(
    vcf_path,
    genes_annotation,
    ref_seqs=None,
    min_freq=0.2,
    min_consensus_freq=0.50,
    curated_models=None,
):
    """Parse VCF file, enforce FILTER, compute FORMAT AD/DP AF, and emit structured coding consequences."""
    variants = []
    if not os.path.exists(vcf_path):
        return variants

    with open(vcf_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            filt = parts[6].strip()
            if filt not in ("PASS", "."):
                continue

            ref_chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            if alt == ".":
                continue
            if "," in alt:
                print(
                    f"Warning: Skipping unsplit multiallelic VCF record at {ref_chrom}:{pos} ({ref}->{alt}) in {vcf_path}",
                    file=sys.stderr,
                )
                continue

            if parts[5] != ".":
                try:
                    qual = round(float(parts[5]), 1)
                except ValueError:
                    qual = None
            else:
                qual = None

            info_str = parts[7]
            info = {}
            for item in info_str.split(";"):
                if not item or item == ".":
                    continue
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                else:
                    info[item] = True

            dp, af = _extract_vcf_dp_and_af(parts, info)
            if dp <= 0 or af < min_freq:
                continue

            if "CONSENSUS" in info:
                is_consensus = True
            elif "SUBCLONAL" in info:
                is_consensus = False
            else:
                is_consensus = af >= min_consensus_freq

            ref_genes = genes_annotation.get(ref_chrom, [])
            raw_matching = [g for g in ref_genes if g["type"] in ["gene", "CDS"] and g["start"] <= pos <= g["end"]]

            # Deduplicate multi-interval CDS parts sharing the same grouping_key/protein so
            # overlapping curated features (e.g. SARS-CoV-2 ORF1ab interval 1 & 2) are annotated once
            seen_features = set()
            matching_genes = []
            for g in raw_matching:
                dedup_key = (
                    g.get("grouping_key") or g.get("protein_id") or g.get("id") or g.get("name"),
                    g.get("strand"),
                )
                if dedup_key not in seen_features:
                    seen_features.add(dedup_key)
                    matching_genes.append(g)

            consequences = []
            if not matching_genes:
                effect_gene = "Intergenic"
                effect_type = "Non-coding"
                aa_change = ""
                status = None
            else:
                effect_type = "Coding"
                gene_labels = []
                resolved_aa = []
                unavailable_statuses = []

                for gene in matching_genes:
                    label = _get_feature_display_name(gene)
                    if label not in gene_labels:
                        gene_labels.append(label)

                    if ref_seqs:
                        res = get_aa_change(
                            pos,
                            ref,
                            alt,
                            gene,
                            ref_seqs.get(ref_chrom, ""),
                            curated_models=curated_models,
                        )
                    else:
                        res = {"aa_change": "", "effect": "unknown", "status": None}

                    consequences.append({
                        "gene": label,
                        "protein_id": gene.get("protein_id") or None,
                        "product": gene.get("product") or None,
                        "aa_change": res["aa_change"],
                        "effect": res.get("effect", "unknown"),
                        "status": res["status"],
                    })

                    if res["aa_change"] and res["aa_change"] != "AA annotation unavailable":
                        if res["aa_change"] not in resolved_aa:
                            resolved_aa.append(res["aa_change"])
                    if res["status"] and res["status"] not in unavailable_statuses:
                        unavailable_statuses.append(res["status"])

                effect_gene = "/".join(gene_labels)
                if resolved_aa:
                    aa_change = "/".join(resolved_aa)
                    status = None
                elif any(c["aa_change"] == "AA annotation unavailable" for c in consequences):
                    aa_change = "AA annotation unavailable"
                    status = unavailable_statuses[0] if unavailable_statuses else None
                else:
                    aa_change = ""
                    status = unavailable_statuses[0] if unavailable_statuses else None

            var_dict = {
                "reference": ref_chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "quality": qual,
                "depth": int(dp),
                "frequency": round(af, 4),
                "is_consensus": bool(is_consensus),
                "gene": effect_gene,
                "type": effect_type,
                "aa_change": aa_change,
                "consequences": consequences,
            }
            if status:
                var_dict["status"] = status
            variants.append(var_dict)

    return variants


def parse_consensus(fasta_path, sample_id=None):
    """Analyze consensus FASTA sequence for all references."""
    results = {}
    if not os.path.exists(fasta_path):
        return results

    consensuses = {}
    current_name = None
    current_seq = []
    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_name:
                    consensuses[current_name] = "".join(current_seq)
                current_name = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_name:
            consensuses[current_name] = "".join(current_seq)

    for name, seq in consensuses.items():
        ref_name = name
        if sample_id and name.startswith(sample_id + "_"):
            ref_name = name[len(sample_id) + 1 :]

        length = len(seq)
        g_count = seq.count("G")
        c_count = seq.count("C")
        n_count = seq.count("N")

        gc_rate = (g_count + c_count) / (length - n_count) if (length - n_count) > 0 else 0.0
        n_rate = n_count / length if length > 0 else 0.0

        results[ref_name] = {
            "length": length,
            "gc_content": round(gc_rate * 100, 2),
            "n_content": round(n_rate * 100, 2),
            "sequence": seq,
        }
    return results


def parse_metadata(metadata_path, processed_sample_ids):
    """Parse sample metadata CSV and validate sample IDs."""
    metadata_map = {}
    if not metadata_path or not os.path.exists(metadata_path):
        return metadata_map

    with open(metadata_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            print("Error: Metadata CSV must contain a 'sample_id' column.", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            s_id = row.get("sample_id")
            if not s_id:
                continue
            metadata_map[s_id] = {k: v for k, v in row.items() if k != "sample_id"}

    processed_set = set(processed_sample_ids)
    metadata_set = set(metadata_map.keys())

    unmatched_metadata = metadata_set - processed_set
    if unmatched_metadata:
        print(
            f"Warning: The following sample IDs in the metadata CSV were not processed in the pipeline run: {', '.join(sorted(unmatched_metadata))}",
            file=sys.stderr,
        )

    unmatched_processed = processed_set - metadata_set
    if unmatched_processed:
        print(
            f"Warning: The following processed sample IDs had no matching metadata entries in the CSV: {', '.join(sorted(unmatched_processed))}",
            file=sys.stderr,
        )

    return metadata_map


_PRIMER_STRICT_RE = re.compile(
    r"[-_](LEFT|RIGHT|FORWARD|REVERSE|FWD|REV|F|R|L)(?:[-_](?:alt)?\d+)?$",
    re.IGNORECASE,
)
_PRE_DIR_ALT_RE = re.compile(r"[-_](?:alt)?\d+$", re.IGNORECASE)


def _primer_side(dir_token, strand):
    if strand == "+":
        return "left"
    if strand == "-":
        return "right"
    if dir_token:
        tok = dir_token.upper()
        if tok in ("LEFT", "FORWARD", "FWD", "F", "L"):
            return "left"
        if tok in ("RIGHT", "REVERSE", "REV", "R"):
            return "right"
    return None


def parse_primer_amplicons(primer_bed_path):
    """Pair primers from BED into trimmed insert amplicons [left_end, right_start) (0-based half-open).

    Supports:
      - ARTIC-style post-direction alt indices: nCoV-2019_1_LEFT_1, nCoV-2019_1_LEFT_alt1, -alt2
      - Pre-direction alt indices: RSV_S1_01_LEFT paired with RSV_S1_RIGHT
      - F/R / FORWARD/REVERSE naming: Pan-CoV-F / Pan-CoV-R
    """
    if not primer_bed_path or not os.path.exists(primer_bed_path):
        return []

    raw_groups = defaultdict(lambda: {"left": [], "right": [], "name": None})

    with open(primer_bed_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            name = parts[3].strip()
            strand = parts[5].strip() if len(parts) >= 6 else "."

            match = _PRIMER_STRICT_RE.search(name)
            dir_token = match.group(1) if match else None
            side = _primer_side(dir_token, strand)
            if side is None:
                continue

            prefix = _PRIMER_STRICT_RE.sub("", name)
            bucket = raw_groups[(chrom, prefix)]
            bucket["name"] = prefix
            bucket[side].append((start, end))

    # Pass 2: Reconcile single-sided alt keys like RSV_S1_01_LEFT -> RSV_S1 with RSV_S1_RIGHT
    merged_groups = defaultdict(lambda: {"left": [], "right": [], "name": None})
    one_sided = []
    for key, val in raw_groups.items():
        if val["left"] and val["right"]:
            merged_groups[key]["left"].extend(val["left"])
            merged_groups[key]["right"].extend(val["right"])
            merged_groups[key]["name"] = key[1]
        else:
            one_sided.append((key, val))

    for (chrom, prefix), val in one_sided:
        alt_prefix = _PRE_DIR_ALT_RE.sub("", prefix)
        target_key = (chrom, prefix)
        if alt_prefix != prefix:
            candidate_key = (chrom, alt_prefix)
            if candidate_key in merged_groups or candidate_key in raw_groups:
                target_key = candidate_key
        merged_groups[target_key]["left"].extend(val["left"])
        merged_groups[target_key]["right"].extend(val["right"])
        merged_groups[target_key]["name"] = target_key[1]

    amplicons = []
    for (chrom, prefix), val in merged_groups.items():
        if not val["left"] or not val["right"]:
            continue
        raw_start = min(p[0] for p in val["left"])
        raw_end = max(p[1] for p in val["right"])
        insert_start = max(p[1] for p in val["left"])
        insert_end = min(p[0] for p in val["right"])
        if insert_start >= insert_end:
            # Degenerate zero/negative insert fallback to full amplicon span
            insert_start, insert_end = raw_start, raw_end

        amplicons.append({
            "chrom": chrom,
            "start": insert_start,
            "end": insert_end,
            "raw_start": raw_start,
            "raw_end": raw_end,
            "name": prefix,
        })

    return amplicons


def parse_explicit_amplicons_bed(amplicons_bed_path):
    """Parse explicit amplicon coordinates BED file."""
    amplicons = []
    if not amplicons_bed_path or not os.path.exists(amplicons_bed_path):
        return amplicons
    with open(amplicons_bed_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                amplicons.append({
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3],
                })
    return amplicons


def _merge_intervals(intervals, max_len):
    """Merge 0-based half-open intervals clamped to [0, max_len)."""
    cleaned = []
    for s, e in intervals:
        cs = max(0, min(max_len, int(s)))
        ce = max(0, min(max_len, int(e)))
        if cs < ce:
            cleaned.append((cs, ce))
    if not cleaned:
        return []
    cleaned.sort()
    merged = [cleaned[0]]
    for cs, ce in cleaned[1:]:
        ps, pe = merged[-1]
        if cs <= pe:
            merged[-1] = (ps, max(pe, ce))
        else:
            merged.append((cs, ce))
    return merged


def calculate_amplicon_coverage(amplicons, depths, min_cov):
    """Calculate per-amplicon coverage stats using and-gated pass criteria."""
    amplicon_metrics = []
    for amp in amplicons:
        chrom = amp["chrom"]
        start = amp["start"]
        end = amp["end"]
        name = amp["name"]

        chrom_depths = depths.get(chrom, [])
        if not chrom_depths:
            continue

        amp_depths = chrom_depths[start : min(end, len(chrom_depths))]
        if not amp_depths:
            continue

        mean_depth = sum(amp_depths) / len(amp_depths)
        min_depth = min(amp_depths)

        sorted_depths = sorted(amp_depths)
        median_depth = sorted_depths[len(sorted_depths) // 2]

        bases_above = sum(1 for d in amp_depths if d >= min_cov)
        pct_above = (bases_above / len(amp_depths)) * 100.0

        status = "pass" if (pct_above >= 90.0 and mean_depth >= min_cov) else "fail"

        amplicon_metrics.append({
            "amplicon_id": name,
            "reference": chrom,
            "start": start + 1,
            "end": min(end, len(chrom_depths)),
            "mean_depth": round(mean_depth, 2),
            "min_depth": min_depth,
            "median_depth": median_depth,
            "percent_bases_above_threshold": round(pct_above, 2),
            "status": status,
        })
    return amplicon_metrics


def parse_and_calculate_amplicons(primer_bed_path, amplicons_bed_path, depths, min_cov):
    """Parse amplicon boundaries (explicit BED preferred for display, else trimmed primer inserts) and compute stats."""
    if amplicons_bed_path and os.path.exists(amplicons_bed_path):
        amplicons = parse_explicit_amplicons_bed(amplicons_bed_path)
        return calculate_amplicon_coverage(amplicons, depths, min_cov)

    if primer_bed_path and os.path.exists(primer_bed_path):
        amplicons = parse_primer_amplicons(primer_bed_path)
        if not amplicons:
            print(
                "Error: Could not reliably pair primers into amplicons from the primer BED file. Please check naming conventions (e.g. using '_LEFT' / '_RIGHT' or '_F' / '_R' suffixes) or supply an explicit amplicon coordinates file using --amplicons-bed.",
                file=sys.stderr,
            )
            sys.exit(1)
        return calculate_amplicon_coverage(amplicons, depths, min_cov)

    return []


def classify_contig_qc(mean_depth, breadth, internal_n_rate, qc_thresholds):
    """Classify per-contig QC tier (PASS/WARN/FAIL) with explicit reason codes."""
    pass_dp = qc_thresholds["pass_depth"]
    pass_br = qc_thresholds["pass_breadth"]
    pass_n = qc_thresholds["pass_max_n"]
    warn_dp = qc_thresholds["warn_depth"]
    warn_br = qc_thresholds["warn_breadth"]
    warn_n = qc_thresholds["warn_max_n"]

    if mean_depth >= pass_dp and breadth >= pass_br and internal_n_rate <= pass_n:
        return "PASS", ["meets_all_pass_thresholds"]

    fail_reasons = []
    if mean_depth < warn_dp:
        fail_reasons.append(f"callable_mean_depth_below_warn ({mean_depth:.2f} < {warn_dp})")
    if breadth < warn_br:
        fail_reasons.append(f"callable_breadth_below_warn ({breadth:.2f}% < {warn_br}%)")
    if internal_n_rate > warn_n:
        fail_reasons.append(f"internal_n_rate_above_warn ({internal_n_rate:.2f}% > {warn_n}%)")
    if fail_reasons:
        return "FAIL", fail_reasons

    warn_reasons = []
    if mean_depth < pass_dp:
        warn_reasons.append(f"callable_mean_depth_below_pass ({mean_depth:.2f} < {pass_dp})")
    if breadth < pass_br:
        warn_reasons.append(f"callable_breadth_below_pass ({breadth:.2f}% < {pass_br}%)")
    if internal_n_rate > pass_n:
        warn_reasons.append(f"internal_n_rate_above_pass ({internal_n_rate:.2f}% > {pass_n}%)")
    return "WARN", warn_reasons


def annotate_callable_targets_and_qc(
    s_data,
    ref_info,
    primer_amplicons,
    explicit_amplicons,
    min_cov,
    mask_depth,
    qc_thresholds,
):
    """Compute per-contig callable-target metrics, internal/terminal Ns, and QC tier classifications."""
    full_depths = s_data.get("full_depths", {})
    coverage = s_data.get("coverage", {})
    consensus = s_data.get("consensus", {})

    primer_by_ref = defaultdict(list)
    for amp in primer_amplicons:
        primer_by_ref[amp["chrom"]].append((amp["start"], amp["end"]))

    explicit_by_ref = defaultdict(list)
    for amp in explicit_amplicons:
        explicit_by_ref[amp["chrom"]].append((amp["start"], amp["end"]))

    best_ref = None
    best_score = (-1.0, -1.0, 1)  # (callable_breadth, callable_mean_depth, -idx)

    for ref_idx, ref in enumerate(ref_info):
        ref_name = ref["name"]
        ref_len = ref["length"]
        ref_seq = ref["sequence"].upper()
        ref_depths = full_depths.get(ref_name, [0] * ref_len)

        if primer_by_ref.get(ref_name):
            qc_basis = "primer_bed"
            target_intervals_0b = _merge_intervals(primer_by_ref[ref_name], ref_len)
        elif explicit_by_ref.get(ref_name):
            qc_basis = "amplicons_bed"
            target_intervals_0b = _merge_intervals(explicit_by_ref[ref_name], ref_len)
        else:
            qc_basis = "whole_contig"
            target_intervals_0b = [(0, ref_len)] if ref_len > 0 else []

        target_mask = [False] * ref_len
        for cs, ce in target_intervals_0b:
            for p in range(cs, ce):
                target_mask[p] = True

        target_depths = [ref_depths[i] for i in range(ref_len) if target_mask[i]]
        callable_len = len(target_depths)
        if callable_len > 0:
            callable_mean_depth = round(sum(target_depths) / callable_len, 2)
            callable_breadth = round((sum(1 for d in target_depths if d >= min_cov) / callable_len) * 100.0, 2)
        else:
            callable_mean_depth = 0.0
            callable_breadth = 0.0

        # Compute masked/N counts on reference coordinates (immune to consensus deletion shifts)
        internal_n_count = 0
        terminal_n_count = 0
        for i in range(ref_len):
            is_n = (ref_depths[i] < mask_depth) or (i < len(ref_seq) and ref_seq[i] == "N")
            if is_n:
                if target_mask[i]:
                    internal_n_count += 1
                else:
                    terminal_n_count += 1

        internal_n_rate = round((internal_n_count / callable_len) * 100.0, 2) if callable_len > 0 else 0.0

        qc_status, qc_reasons = classify_contig_qc(
            callable_mean_depth,
            callable_breadth,
            internal_n_rate,
            qc_thresholds,
        )

        if ref_name in coverage:
            coverage[ref_name]["callable_target_length"] = callable_len
            coverage[ref_name]["callable_target_intervals"] = [[cs + 1, ce] for cs, ce in target_intervals_0b]
            coverage[ref_name]["callable_mean_depth"] = callable_mean_depth
            coverage[ref_name]["callable_breadth_min_cov"] = callable_breadth
            coverage[ref_name]["qc_basis"] = qc_basis
            coverage[ref_name]["qc_status"] = qc_status
            coverage[ref_name]["qc_reasons"] = qc_reasons

        if ref_name in consensus:
            consensus[ref_name]["internal_n_content"] = internal_n_rate
            consensus[ref_name]["internal_n_count"] = internal_n_count
            consensus[ref_name]["terminal_n_count"] = terminal_n_count

        score = (callable_breadth, callable_mean_depth, -ref_idx)
        if score > best_score:
            best_score = score
            best_ref = ref_name

    if best_ref and best_ref in coverage:
        s_data["primary_reference"] = best_ref
        s_data["qc_status"] = coverage[best_ref]["qc_status"]
        s_data["qc_reasons"] = coverage[best_ref]["qc_reasons"]
    else:
        s_data["primary_reference"] = ref_info[0]["name"] if ref_info else None
        s_data["qc_status"] = "FAIL"
        s_data["qc_reasons"] = ["no_reference_coverage"]


def calculate_distance_matrices(samples_data, ref_info, min_cov, min_consensus_freq=0.50):
    """Calculate coordinate-aligned SNP-only pairwise distances, raw SNP counts, and compared_sites.

    Assigns None (JSON null) instead of 1.0 when zero mutually covered sites exist.
    """
    distance_matrices = {}
    snp_distance_matrices = {}
    compared_sites_matrices = {}
    distance_matrices_by_name = {}

    samples_list = list(samples_data.values())
    num_samples = len(samples_list)
    distance_matrix_order = [s["sample_id"] for s in samples_list]

    for ref in ref_info:
        ref_name = ref["name"]
        ref_len = ref["length"]
        ref_seq = ref["sequence"]

        matrix = [[0.0] * num_samples for _ in range(num_samples)]
        snp_matrix = [[0] * num_samples for _ in range(num_samples)]
        sites_matrix = [[0] * num_samples for _ in range(num_samples)]

        per_sample_vars = []
        per_sample_dels = []
        per_sample_depths = []

        for s in samples_list:
            depths_s = s.get("full_depths", {}).get(ref_name, [0] * ref_len)
            vars_s = {
                v["position"]: v["alt"]
                for v in s.get("variants", [])
                if v["reference"] == ref_name
                and len(v["ref"]) == 1
                and len(v["alt"]) == 1
                and v.get("frequency", 1.0) >= min_consensus_freq
            }
            del_positions_s = set()
            for v in s.get("variants", []):
                if v["reference"] == ref_name and len(v["ref"]) > len(v["alt"]) and v.get("frequency", 1.0) >= min_consensus_freq:
                    start_del = v["position"] + len(v["alt"])
                    end_del = v["position"] + len(v["ref"]) - 1
                    del_positions_s.update(range(start_del, end_del + 1))

            per_sample_depths.append(depths_s)
            per_sample_vars.append(vars_s)
            per_sample_dels.append(del_positions_s)

        for i in range(num_samples):
            depths_A = per_sample_depths[i]
            vars_A = per_sample_vars[i]
            del_A = per_sample_dels[i]

            self_sites = sum(
                1
                for pos_idx in range(ref_len)
                if (pos_idx + 1) not in del_A and depths_A[pos_idx] >= min_cov
            )
            sites_matrix[i][i] = self_sites

            for j in range(i + 1, num_samples):
                depths_B = per_sample_depths[j]
                vars_B = per_sample_vars[j]
                del_B = per_sample_dels[j]

                diffs = 0
                compared = 0

                for pos_idx in range(ref_len):
                    pos = pos_idx + 1
                    if pos in del_A or pos in del_B:
                        continue
                    if depths_A[pos_idx] >= min_cov and depths_B[pos_idx] >= min_cov:
                        compared += 1
                        base_A = vars_A.get(pos, ref_seq[pos_idx])
                        base_B = vars_B.get(pos, ref_seq[pos_idx])
                        if base_A.upper() != base_B.upper():
                            diffs += 1

                if compared > 0:
                    dist = diffs / compared
                    snp_val = diffs
                else:
                    dist = None
                    snp_val = None

                matrix[i][j] = dist
                matrix[j][i] = dist
                snp_matrix[i][j] = snp_val
                snp_matrix[j][i] = snp_val
                sites_matrix[i][j] = compared
                sites_matrix[j][i] = compared

        distance_matrices[ref_name] = matrix
        snp_distance_matrices[ref_name] = snp_matrix
        compared_sites_matrices[ref_name] = sites_matrix

        ref_by_name = {}
        for i, sA_id in enumerate(distance_matrix_order):
            ref_by_name[sA_id] = {}
            for j, sB_id in enumerate(distance_matrix_order):
                ref_by_name[sA_id][sB_id] = matrix[i][j]
        distance_matrices_by_name[ref_name] = ref_by_name

    return (
        distance_matrices,
        distance_matrix_order,
        distance_matrices_by_name,
        snp_distance_matrices,
        compared_sites_matrices,
    )


def parse_versions_files(versions_arg):
    """Parse tool: version pairs from process versions_*.yml files."""
    versions = {"python": sys.version.split()[0]}
    if versions_arg:
        for path in [p.strip() for p in versions_arg.split(",") if p.strip()]:
            if not os.path.exists(path):
                continue
            with open(path, "r") as f:
                for line in f:
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#") or ":" not in stripped:
                        continue
                    k, v = stripped.split(":", 1)
                    k = k.strip()
                    v = v.strip().strip("'\"")
                    if k and v:
                        versions[k] = v
    if "samtools" not in versions:
        try:
            res = subprocess.run(["samtools", "--version"], capture_output=True, text=True, check=True)
            first = res.stdout.splitlines()[0] if res.stdout else ""
            if first.startswith("samtools "):
                versions["samtools"] = first.split()[1]
        except Exception:
            pass
    return dict(sorted(versions.items()))


def load_run_config(args, min_consensus_depth):
    """Load run_config JSON staged by Nextflow or synthesize from CLI arguments."""
    default_cfg = {
        "run_id": args.run_id,
        "platform": "unknown",
        "min_coverage": args.min_coverage,
        "min_freq": args.min_freq,
        "min_consensus_depth": min_consensus_depth,
        "min_consensus_freq": args.min_consensus_freq,
        "min_map_qual": args.min_map_qual,
        "min_base_qual": args.min_base_qual,
        "min_variant_qual": None,
        "reference": os.path.basename(args.reference) if args.reference else None,
        "gff": os.path.basename(args.gff) if args.gff else None,
        "primer_bed": os.path.basename(args.primers) if args.primers else None,
        "amplicons_bed": os.path.basename(args.amplicons_bed) if args.amplicons_bed else None,
        "metadata": os.path.basename(args.metadata) if args.metadata else None,
        "cds_models": os.path.basename(args.cds_models) if args.cds_models else None,
        "qc_thresholds": {
            "pass_depth": args.qc_pass_depth,
            "pass_breadth": args.qc_pass_breadth,
            "pass_max_n": args.qc_pass_max_n,
            "warn_depth": args.qc_warn_depth,
            "warn_breadth": args.qc_warn_breadth,
            "warn_max_n": args.qc_warn_max_n,
        },
    }
    if args.run_config and os.path.exists(args.run_config):
        with open(args.run_config, "r") as f:
            loaded = json.load(f)
        if isinstance(loaded, dict):
            default_cfg.update(loaded)
    for int_k in (
        "min_coverage",
        "min_consensus_depth",
        "min_map_qual",
        "min_base_qual",
        "min_variant_qual",
        "primer_tolerance",
        "min_clipped_len",
        "max_mpileup_depth",
        "downsample_seed",
        "ont_min_mean_qual",
        "ont_min_length",
        "ont_max_length",
    ):
        if default_cfg.get(int_k) is not None:
            default_cfg[int_k] = int(default_cfg[int_k])
    for flt_k in ("min_freq", "min_consensus_freq", "min_iupac_freq", "max_iupac_freq"):
        if default_cfg.get(flt_k) is not None:
            default_cfg[flt_k] = float(default_cfg[flt_k])
    if isinstance(default_cfg.get("qc_thresholds"), dict):
        qt = default_cfg["qc_thresholds"]
        for k in ("pass_depth", "warn_depth"):
            if qt.get(k) is not None:
                qt[k] = int(qt[k])
        for k in ("pass_breadth", "pass_max_n", "warn_breadth", "warn_max_n"):
            if qt.get(k) is not None:
                qt[k] = float(qt[k])
    return default_cfg


def _is_num(val):
    return isinstance(val, (int, float)) and not isinstance(val, bool)


def validate_report_schema(report):
    """
    Zero-dependency validator checking report against the Milestone 2
    schema/viromics_run_report.schema.json contract.
    """
    if not isinstance(report, dict):
        raise ValueError("Schema Error: Report must be a JSON object dictionary.")

    required_top = [
        "schema_version",
        "run_id",
        "run_config",
        "software_versions",
        "min_coverage_threshold",
        "references",
        "genes",
        "distance_matrices",
        "distance_matrix_order",
        "distance_matrices_by_name",
        "snp_distance_matrices",
        "compared_sites_matrices",
        "samples",
    ]
    for key in required_top:
        if key not in report:
            raise ValueError(f"Schema Error: Report missing required top-level key '{key}'.")

    if not isinstance(report["schema_version"], str) or not report["schema_version"]:
        raise ValueError("Schema Error: 'schema_version' must be a non-empty string.")
    if not isinstance(report["run_id"], str) or not report["run_id"]:
        raise ValueError("Schema Error: 'run_id' must be a non-empty string.")
    if not isinstance(report["run_config"], dict):
        raise ValueError("Schema Error: 'run_config' must be an object dictionary.")
    for cfg_k in ("platform", "min_coverage", "min_freq", "min_consensus_depth", "min_consensus_freq", "min_map_qual", "min_base_qual", "qc_thresholds"):
        if cfg_k not in report["run_config"]:
            raise ValueError(f"Schema Error: 'run_config' missing required key '{cfg_k}'.")
    qc_t = report["run_config"]["qc_thresholds"]
    if not isinstance(qc_t, dict):
        raise ValueError("Schema Error: 'run_config.qc_thresholds' must be a dictionary.")
    for tk in ("pass_depth", "pass_breadth", "pass_max_n", "warn_depth", "warn_breadth", "warn_max_n"):
        if tk not in qc_t or not _is_num(qc_t[tk]):
            raise ValueError(f"Schema Error: 'run_config.qc_thresholds' missing numeric key '{tk}'.")

    if not isinstance(report["software_versions"], dict):
        raise ValueError("Schema Error: 'software_versions' must be a dictionary.")
    for k, v in report["software_versions"].items():
        if not isinstance(k, str) or not isinstance(v, str):
            raise ValueError(f"Schema Error: 'software_versions' entry '{k}' must map string->string.")

    if not isinstance(report["min_coverage_threshold"], int) or isinstance(report["min_coverage_threshold"], bool):
        raise ValueError("Schema Error: 'min_coverage_threshold' must be an integer.")
    if not isinstance(report["references"], list):
        raise ValueError("Schema Error: 'references' must be a list.")
    if not isinstance(report["genes"], dict):
        raise ValueError("Schema Error: 'genes' must be an object dictionary.")
    if not isinstance(report["distance_matrix_order"], list):
        raise ValueError("Schema Error: 'distance_matrix_order' must be a list.")
    for idx, val in enumerate(report["distance_matrix_order"]):
        if not isinstance(val, str):
            raise ValueError(f"Schema Error: 'distance_matrix_order' entry at index {idx} must be a string.")

    for mat_key in ("distance_matrices", "snp_distance_matrices", "compared_sites_matrices"):
        if not isinstance(report[mat_key], dict):
            raise ValueError(f"Schema Error: '{mat_key}' must be a dictionary.")
        for ref_name, rows in report[mat_key].items():
            if not isinstance(rows, list):
                raise ValueError(f"Schema Error: '{mat_key}[{ref_name}]' must be a 2D list.")
            for r in rows:
                if not isinstance(r, list):
                    raise ValueError(f"Schema Error: '{mat_key}[{ref_name}]' row must be a list.")
                for val in r:
                    if mat_key == "compared_sites_matrices":
                        if not isinstance(val, int) or isinstance(val, bool) or val < 0:
                            raise ValueError(f"Schema Error: '{mat_key}[{ref_name}]' entry must be a non-negative integer.")
                    elif mat_key == "snp_distance_matrices":
                        if val is not None and (not isinstance(val, int) or isinstance(val, bool) or val < 0):
                            raise ValueError(f"Schema Error: '{mat_key}[{ref_name}]' entry must be integer or null.")
                    else:
                        if val is not None and (not _is_num(val) or val < 0):
                            raise ValueError(f"Schema Error: 'distance_matrices[{ref_name}]' entry must be number or null.")

    if not isinstance(report["distance_matrices_by_name"], dict):
        raise ValueError("Schema Error: 'distance_matrices_by_name' must be an object dictionary.")
    for ref_name, matrix_dict in report["distance_matrices_by_name"].items():
        if not isinstance(matrix_dict, dict):
            raise ValueError(f"Schema Error: 'distance_matrices_by_name' for ref '{ref_name}' must be a dictionary.")
        for sA_id, sB_dict in matrix_dict.items():
            if not isinstance(sB_dict, dict):
                raise ValueError(f"Schema Error: 'distance_matrices_by_name' entry for ref '{ref_name}', sample '{sA_id}' must be a dictionary.")
            for sB_id, dist in sB_dict.items():
                if dist is not None and not _is_num(dist):
                    raise ValueError(f"Schema Error: 'distance_matrices_by_name' distance for ref '{ref_name}', samples '{sA_id}' vs '{sB_id}' must be a number or null.")

    if not isinstance(report["samples"], list):
        raise ValueError("Schema Error: 'samples' must be a list.")

    for idx, ref in enumerate(report["references"]):
        if not isinstance(ref, dict):
            raise ValueError(f"Schema Error: Reference entry at index {idx} is not a dictionary.")
        for key in ["name", "length", "sequence"]:
            if key not in ref:
                raise ValueError(f"Schema Error: Reference at index {idx} missing key '{key}'.")
        if not isinstance(ref["name"], str):
            raise ValueError(f"Schema Error: Reference at index {idx} 'name' must be a string.")
        if not isinstance(ref["length"], int) or isinstance(ref["length"], bool):
            raise ValueError(f"Schema Error: Reference at index {idx} 'length' must be an integer.")
        if not isinstance(ref["sequence"], str):
            raise ValueError(f"Schema Error: Reference at index {idx} 'sequence' must be a string.")

    valid_tiers = {"PASS", "WARN", "FAIL"}
    for idx, s in enumerate(report["samples"]):
        if not isinstance(s, dict):
            raise ValueError(f"Schema Error: Sample at index {idx} is not a dictionary.")
        required_sample = ["sample_id", "qc_status", "qc_reasons", "primary_reference", "qc", "coverage", "variants", "consensus", "metadata", "amplicons"]
        for key in required_sample:
            if key not in s:
                raise ValueError(f"Schema Error: Sample at index {idx} (ID: {s.get('sample_id', 'unknown')}) missing required key '{key}'.")

        s_id = s["sample_id"]
        if not isinstance(s_id, str):
            raise ValueError(f"Schema Error: Sample at index {idx} 'sample_id' must be a string.")
        if s["qc_status"] not in valid_tiers:
            raise ValueError(f"Schema Error: Sample '{s_id}' invalid 'qc_status': {s['qc_status']}.")
        if not isinstance(s["qc_reasons"], list):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'qc_reasons' must be a list.")
        if s["primary_reference"] is not None and not isinstance(s["primary_reference"], str):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'primary_reference' must be string or null.")

        qc = s["qc"]
        if not isinstance(qc, dict):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'qc' must be a dictionary.")
        for key in ["raw_reads", "raw_bases", "trimmed_reads", "trimmed_bases"]:
            if key not in qc or not isinstance(qc[key], int) or isinstance(qc[key], bool):
                raise ValueError(f"Schema Error: Sample '{s_id}' QC key '{key}' must be an integer.")
        for key in ["q30_rate_raw", "q30_rate_trimmed", "duplication_rate"]:
            if key not in qc:
                raise ValueError(f"Schema Error: Sample '{s_id}' QC missing key '{key}'.")
            if qc[key] is not None and not _is_num(qc[key]):
                raise ValueError(f"Schema Error: Sample '{s_id}' QC key '{key}' must be number or null.")

        coverage = s["coverage"]
        if not isinstance(coverage, dict):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'coverage' must be a dictionary.")
        for ref_name, cov_data in coverage.items():
            if not isinstance(cov_data, dict):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage data for ref '{ref_name}' must be a dictionary.")
            required_cov = [
                "mean_depth",
                "coverage_breadth_10x",
                "coverage_breadth_30x",
                "coverage_breadth_min_cov",
                "mapped_reads",
                "coverage_profile",
                "low_coverage_intervals",
                "callable_target_length",
                "callable_target_intervals",
                "callable_mean_depth",
                "callable_breadth_min_cov",
                "qc_basis",
                "qc_status",
                "qc_reasons",
            ]
            for key in required_cov:
                if key not in cov_data:
                    raise ValueError(f"Schema Error: Sample '{s_id}' coverage for ref '{ref_name}' missing key '{key}'.")
            if cov_data["qc_status"] not in valid_tiers:
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage ref '{ref_name}' invalid qc_status.")
            if not isinstance(cov_data["low_coverage_intervals"], list):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage ref '{ref_name}' 'low_coverage_intervals' must be a list.")

        for v_idx, var in enumerate(s["variants"]):
            if not isinstance(var, dict):
                raise ValueError(f"Schema Error: Sample '{s_id}' variant {v_idx} must be a dict.")
            for vkey in ["reference", "position", "ref", "alt", "quality", "depth", "frequency", "is_consensus", "gene", "type", "aa_change", "consequences"]:
                if vkey not in var:
                    raise ValueError(f"Schema Error: Sample '{s_id}' variant {v_idx} missing '{vkey}'.")
            if var["quality"] is not None and not _is_num(var["quality"]):
                raise ValueError(f"Schema Error: Sample '{s_id}' variant {v_idx} 'quality' must be number or null.")
            if not isinstance(var["consequences"], list):
                raise ValueError(f"Schema Error: Sample '{s_id}' variant {v_idx} 'consequences' must be a list.")

        consensus = s["consensus"]
        if not isinstance(consensus, dict):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'consensus' must be a dictionary.")
        for ref_name, cons_data in consensus.items():
            if not isinstance(cons_data, dict):
                raise ValueError(f"Schema Error: Sample '{s_id}' consensus data for ref '{ref_name}' must be a dictionary.")
            required_cons = ["length", "gc_content", "n_content", "internal_n_content", "internal_n_count", "terminal_n_count", "sequence"]
            for key in required_cons:
                if key not in cons_data:
                    raise ValueError(f"Schema Error: Sample '{s_id}' consensus for ref '{ref_name}' missing key '{key}'.")


def main():
    args = parse_args()
    min_consensus_depth = args.min_consensus_depth if args.min_consensus_depth is not None else args.min_coverage

    ref_info = parse_reference(args.reference)
    ref_lengths = {ref["name"]: ref["length"] for ref in ref_info}
    ref_seqs = {ref["name"]: ref["sequence"] for ref in ref_info}

    genes = parse_gff(args.gff)
    curated_models = load_curated_models_cached(args.cds_models)

    fastp_files = sorted(args.fastp.split(","))
    bam_files = sorted(args.bams.split(","))
    vcf_files = sorted(args.vcfs.split(","))
    consensus_files = sorted(args.consensuses.split(","))

    samples_data = {}
    for fastp_f in fastp_files:
        sample_id = os.path.basename(fastp_f).replace("_fastp.json", "")
        samples_data[sample_id] = {
            "sample_id": sample_id,
            "qc": parse_fastp(fastp_f),
        }

    for bam_f in bam_files:
        sample_id = os.path.basename(bam_f).replace("_trimmed.bam", "").replace("_raw.bam", "")
        if sample_id in samples_data:
            cov, full_depths = calculate_coverage(
                bam_f,
                ref_lengths,
                min_cov=args.min_coverage,
                min_consensus_depth=min_consensus_depth,
                min_map_qual=args.min_map_qual,
                min_base_qual=args.min_base_qual,
            )
            samples_data[sample_id]["coverage"] = cov
            samples_data[sample_id]["full_depths"] = full_depths

    for vcf_f in vcf_files:
        sample_id = os.path.basename(vcf_f).replace(".vcf", "")
        if sample_id in samples_data:
            variants = parse_vcf(
                vcf_f,
                genes,
                ref_seqs,
                min_freq=args.min_freq,
                min_consensus_freq=args.min_consensus_freq,
                curated_models=curated_models,
            )
            samples_data[sample_id]["variants"] = variants

    for cons_f in consensus_files:
        sample_id = os.path.basename(cons_f).replace("_consensus.fasta", "")
        if sample_id in samples_data:
            samples_data[sample_id]["consensus"] = parse_consensus(cons_f, sample_id)

    metadata_map = parse_metadata(args.metadata, list(samples_data.keys()))

    primer_amplicons = parse_primer_amplicons(args.primers) if args.primers else []
    explicit_amplicons = parse_explicit_amplicons_bed(args.amplicons_bed) if args.amplicons_bed else []
    qc_thresholds = {
        "pass_depth": args.qc_pass_depth,
        "pass_breadth": args.qc_pass_breadth,
        "pass_max_n": args.qc_pass_max_n,
        "warn_depth": args.qc_warn_depth,
        "warn_breadth": args.qc_warn_breadth,
        "warn_max_n": args.qc_warn_max_n,
    }

    for sample_id, s_data in samples_data.items():
        s_data.setdefault("coverage", {})
        s_data.setdefault("variants", [])
        s_data.setdefault("consensus", {})
        s_data["metadata"] = metadata_map.get(sample_id, {})
        s_data["amplicons"] = []
        if (args.primers or args.amplicons_bed) and "full_depths" in s_data:
            s_data["amplicons"] = parse_and_calculate_amplicons(
                args.primers, args.amplicons_bed, s_data["full_depths"], args.min_coverage
            )
        annotate_callable_targets_and_qc(
            s_data,
            ref_info,
            primer_amplicons,
            explicit_amplicons,
            min_cov=args.min_coverage,
            mask_depth=min_consensus_depth,
            qc_thresholds=qc_thresholds,
        )

    (
        distance_matrices,
        distance_matrix_order,
        distance_matrices_by_name,
        snp_distance_matrices,
        compared_sites_matrices,
    ) = calculate_distance_matrices(
        samples_data,
        ref_info,
        min_cov=min_consensus_depth,
        min_consensus_freq=args.min_consensus_freq,
    )

    for s in samples_data.values():
        if "full_depths" in s:
            del s["full_depths"]

    run_config = load_run_config(args, min_consensus_depth)
    software_versions = parse_versions_files(args.versions)

    report = {
        "schema_version": SCHEMA_VERSION,
        "run_id": args.run_id,
        "run_config": run_config,
        "software_versions": software_versions,
        "min_coverage_threshold": args.min_coverage,
        "references": ref_info,
        "genes": genes,
        "distance_matrices": distance_matrices,
        "distance_matrix_order": distance_matrix_order,
        "distance_matrices_by_name": distance_matrices_by_name,
        "snp_distance_matrices": snp_distance_matrices,
        "compared_sites_matrices": compared_sites_matrices,
        "samples": list(samples_data.values()),
    }

    validate_report_schema(report)

    with open(args.output, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Successfully generated Viromics Run Report: {args.output}")


if __name__ == "__main__":
    main()
