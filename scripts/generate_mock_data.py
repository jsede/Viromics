#!/usr/bin/env python3
"""Generate synthetic multi-strain RSV-A/B Viromics run report conforming to schema v2.0.0."""

import json
import math
import os
import random
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "pipeline", "bin"))
import parse_metrics  # noqa: E402


def generate_amplicon_profile(length, num_bins=500, mean_depth=300, dropout=False, dropout_regions=None):
    """Simulates realistic amplicon coverage depth with peaks and troughs."""
    profile = []
    bin_size = length / num_bins

    num_amplicons = 20
    amp_positions = [i * (length / num_amplicons) for i in range(num_amplicons)]
    amp_heights = [random.randint(int(mean_depth * 0.5), int(mean_depth * 1.5)) for _ in range(num_amplicons)]

    if dropout and dropout_regions:
        for idx in range(num_amplicons):
            pos = amp_positions[idx]
            for start, end in dropout_regions:
                if start <= pos <= end:
                    amp_heights[idx] = random.randint(0, 3)

    for i in range(num_bins):
        pos = i * bin_size
        depth = 2
        for amp_pos, amp_height in zip(amp_positions, amp_heights):
            dist = pos - amp_pos
            width = length / 30
            weight = math.exp(-(dist ** 2) / (2 * (width ** 2)))
            depth += amp_height * weight

        depth = max(1, depth * random.uniform(0.85, 1.15))
        if dropout and dropout_regions:
            for start, end in dropout_regions:
                if start <= pos <= end:
                    depth = max(0, depth * random.uniform(0.01, 0.05))

        profile.append({
            "pos": int(pos) + 1,
            "depth": round(depth),
        })
    return profile


def generate_consensus_sequence(ref_seq, variants, ref_name, dropout=False, dropout_regions=None):
    seq_list = list(ref_seq)
    for var in variants:
        if var["reference"] == ref_name and var.get("is_consensus", True):
            pos_idx = var["position"] - 1
            if 0 <= pos_idx < len(seq_list):
                seq_list[pos_idx] = var["alt"]
    if dropout and dropout_regions:
        for start, end in dropout_regions:
            for i in range(start - 1, end):
                if 0 <= i < len(seq_list):
                    seq_list[i] = "N"
    return "".join(seq_list)


def _make_variant(ref_name, pos, ref, alt, qual, depth, freq, gene_name, product, ref_seq, gene_dict):
    aa_res = parse_metrics.get_aa_change(pos, ref, alt, gene_dict, ref_seq)
    return {
        "reference": ref_name,
        "position": pos,
        "ref": ref,
        "alt": alt,
        "quality": qual,
        "depth": depth,
        "frequency": freq,
        "is_consensus": freq >= 0.5,
        "gene": gene_name,
        "type": "Coding",
        "aa_change": aa_res["aa_change"],
        "consequences": [
            {
                "gene": gene_name,
                "protein_id": None,
                "product": product,
                "aa_change": aa_res["aa_change"],
                "effect": aa_res.get("effect", "unknown"),
                "status": aa_res["status"],
            }
        ],
    }


def _coverage_entry(mean_depth, b10, b30, bmin, mapped_reads, profile, low_cov_intervals, ref_len, qc_thresholds):
    qc_status, qc_reasons = parse_metrics.classify_contig_qc(mean_depth, bmin, round((sum(e - s + 1 for s, e in low_cov_intervals) / ref_len) * 100.0, 2), qc_thresholds)
    return {
        "mean_depth": mean_depth,
        "coverage_breadth_10x": b10,
        "coverage_breadth_30x": b30,
        "coverage_breadth_min_cov": bmin,
        "mapped_reads": mapped_reads,
        "coverage_profile": profile,
        "low_coverage_intervals": low_cov_intervals,
        "callable_target_length": ref_len,
        "callable_target_intervals": [[1, ref_len]],
        "callable_mean_depth": mean_depth,
        "callable_breadth_min_cov": bmin,
        "qc_basis": "whole_contig",
        "qc_status": qc_status,
        "qc_reasons": qc_reasons,
    }


def _consensus_entry(ref_len, gc_pct, cons_seq, low_cov_intervals):
    n_count = cons_seq.count("N")
    low_cov_bases = sum(e - s + 1 for s, e in low_cov_intervals)
    int_n = max(n_count, low_cov_bases) if low_cov_bases == ref_len else n_count
    n_pct = round((n_count / ref_len) * 100.0, 2)
    int_n_pct = round((int_n / ref_len) * 100.0, 2)
    return {
        "length": ref_len,
        "gc_content": gc_pct,
        "n_content": n_pct,
        "internal_n_content": int_n_pct,
        "internal_n_count": int_n,
        "terminal_n_count": 0,
        "sequence": cons_seq,
    }


def main():
    print("Generating Multi-Strain RSV-A/B Mock Run Report (schema v2.0.0)...")
    random.seed(42)

    ref_length = 15200
    bases = ["A", "C", "G", "T"]
    ref_seq_a = "".join(random.choices(bases, k=ref_length))
    ref_seq_b = "".join(random.choices(bases, k=ref_length))

    references = [
        {"name": "RSV-A_strain_Wuhan_2026", "length": ref_length, "sequence": ref_seq_a},
        {"name": "RSV-B_strain_Melbourne_2026", "length": ref_length, "sequence": ref_seq_b},
    ]

    genes = {
        "RSV-A_strain_Wuhan_2026": [
            {"name": "N_Protein", "type": "gene", "start": 100, "end": 1200, "strand": "+", "phase": 0, "product": "Nucleoprotein"},
            {"name": "F_Protein", "type": "gene", "start": 1500, "end": 3500, "strand": "+", "phase": 0, "product": "Fusion glycoprotein"},
            {"name": "G_Protein", "type": "gene", "start": 4000, "end": 6000, "strand": "+", "phase": 0, "product": "Attachment glycoprotein"},
            {"name": "L_Protein", "type": "gene", "start": 7000, "end": 14500, "strand": "+", "phase": 0, "product": "RNA-directed RNA polymerase"},
        ],
        "RSV-B_strain_Melbourne_2026": [
            {"name": "N_Protein", "type": "gene", "start": 100, "end": 1200, "strand": "+", "phase": 0, "product": "Nucleoprotein"},
            {"name": "F_Protein", "type": "gene", "start": 1500, "end": 3500, "strand": "+", "phase": 0, "product": "Fusion glycoprotein"},
            {"name": "G_Protein", "type": "gene", "start": 3800, "end": 5800, "strand": "+", "phase": 0, "product": "Attachment glycoprotein"},
            {"name": "L_Protein", "type": "gene", "start": 7000, "end": 14500, "strand": "+", "phase": 0, "product": "RNA-directed RNA polymerase"},
        ],
    }

    qc_thresholds = {
        "pass_depth": 100,
        "pass_breadth": 95.0,
        "pass_max_n": 5.0,
        "warn_depth": 25,
        "warn_breadth": 80.0,
        "warn_max_n": 10.0,
    }

    g_a = {g["name"]: g for g in genes["RSV-A_strain_Wuhan_2026"]}
    g_b = {g["name"]: g for g in genes["RSV-B_strain_Melbourne_2026"]}

    def _alt_different(ref_b):
        return "T" if ref_b != "T" else "C"

    samples = []

    # Sample 1: Pure RSV-A Infection (PASS on RSV-A, FAIL on RSV-B)
    s1_vars = [
        _make_variant("RSV-A_strain_Wuhan_2026", 520, ref_seq_a[519], _alt_different(ref_seq_a[519]), 99.0, 310, 0.99, "N_Protein", "Nucleoprotein", ref_seq_a, g_a["N_Protein"]),
        _make_variant("RSV-A_strain_Wuhan_2026", 2105, ref_seq_a[2104], _alt_different(ref_seq_a[2104]), 95.0, 280, 0.97, "F_Protein", "Fusion glycoprotein", ref_seq_a, g_a["F_Protein"]),
        _make_variant("RSV-A_strain_Wuhan_2026", 4820, ref_seq_a[4819], _alt_different(ref_seq_a[4819]), 88.0, 410, 0.94, "G_Protein", "Attachment glycoprotein", ref_seq_a, g_a["G_Protein"]),
    ]
    s1_cons_a = generate_consensus_sequence(ref_seq_a, s1_vars, "RSV-A_strain_Wuhan_2026")
    s1_cons_b = generate_consensus_sequence(ref_seq_b, s1_vars, "RSV-B_strain_Melbourne_2026", dropout=True, dropout_regions=[(1, ref_length)])
    s1_cov_a_entry = _coverage_entry(352.4, 99.8, 99.2, 99.8, 105700, generate_amplicon_profile(ref_length, mean_depth=350), [], ref_length, qc_thresholds)
    s1_cov_b_entry = _coverage_entry(0.8, 0.0, 0.0, 0.0, 240, generate_amplicon_profile(ref_length, mean_depth=1), [[1, ref_length]], ref_length, qc_thresholds)
    samples.append({
        "sample_id": "RSV-S001_A_Positive",
        "qc_status": s1_cov_a_entry["qc_status"],
        "qc_reasons": s1_cov_a_entry["qc_reasons"],
        "primary_reference": "RSV-A_strain_Wuhan_2026",
        "qc": {
            "raw_reads": 128042,
            "raw_bases": 38412600,
            "trimmed_reads": 112040,
            "trimmed_bases": 33612000,
            "q30_rate_raw": 0.948,
            "q30_rate_trimmed": 0.969,
            "duplication_rate": 0.0642,
            "mean_read_length_raw": 150.0,
            "mean_read_length_trimmed": 148.5,
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": s1_cov_a_entry,
            "RSV-B_strain_Melbourne_2026": s1_cov_b_entry,
        },
        "variants": s1_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": _consensus_entry(ref_length, 48.2, s1_cons_a, []),
            "RSV-B_strain_Melbourne_2026": _consensus_entry(ref_length, 45.1, s1_cons_b, [[1, ref_length]]),
        },
        "metadata": {"host": "Human", "collection_date": "2026-05-20", "location": "Sydney"},
        "amplicons": [],
    })

    # Sample 2: Pure RSV-B Infection (PASS on RSV-B, FAIL on RSV-A)
    s2_vars = [
        _make_variant("RSV-B_strain_Melbourne_2026", 612, ref_seq_b[611], _alt_different(ref_seq_b[611]), 99.0, 240, 0.99, "N_Protein", "Nucleoprotein", ref_seq_b, g_b["N_Protein"]),
        _make_variant("RSV-B_strain_Melbourne_2026", 4105, ref_seq_b[4104], _alt_different(ref_seq_b[4104]), 92.0, 320, 0.95, "G_Protein", "Attachment glycoprotein", ref_seq_b, g_b["G_Protein"]),
    ]
    s2_cons_a = generate_consensus_sequence(ref_seq_a, s2_vars, "RSV-A_strain_Wuhan_2026", dropout=True, dropout_regions=[(1, ref_length)])
    s2_cons_b = generate_consensus_sequence(ref_seq_b, s2_vars, "RSV-B_strain_Melbourne_2026")
    s2_cov_a_entry = _coverage_entry(0.6, 0.0, 0.0, 0.0, 180, generate_amplicon_profile(ref_length, mean_depth=1), [[1, ref_length]], ref_length, qc_thresholds)
    s2_cov_b_entry = _coverage_entry(281.2, 99.5, 98.4, 99.5, 84300, generate_amplicon_profile(ref_length, mean_depth=280), [], ref_length, qc_thresholds)
    samples.append({
        "sample_id": "RSV-S002_B_Positive",
        "qc_status": s2_cov_b_entry["qc_status"],
        "qc_reasons": s2_cov_b_entry["qc_reasons"],
        "primary_reference": "RSV-B_strain_Melbourne_2026",
        "qc": {
            "raw_reads": 98402,
            "raw_bases": 29520600,
            "trimmed_reads": 86240,
            "trimmed_bases": 25872000,
            "q30_rate_raw": 0.939,
            "q30_rate_trimmed": 0.961,
            "duplication_rate": 0.0821,
            "mean_read_length_raw": 150.0,
            "mean_read_length_trimmed": 147.8,
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": s2_cov_a_entry,
            "RSV-B_strain_Melbourne_2026": s2_cov_b_entry,
        },
        "variants": s2_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": _consensus_entry(ref_length, 48.2, s2_cons_a, [[1, ref_length]]),
            "RSV-B_strain_Melbourne_2026": _consensus_entry(ref_length, 44.8, s2_cons_b, []),
        },
        "metadata": {"host": "Human", "collection_date": "2026-05-21", "location": "Melbourne"},
        "amplicons": [],
    })

    # Sample 3: Co-infection / Co-detection (PASS on both!)
    s3_vars = [
        _make_variant("RSV-A_strain_Wuhan_2026", 520, ref_seq_a[519], _alt_different(ref_seq_a[519]), 99.0, 180, 0.99, "N_Protein", "Nucleoprotein", ref_seq_a, g_a["N_Protein"]),
        _make_variant("RSV-B_strain_Melbourne_2026", 612, ref_seq_b[611], _alt_different(ref_seq_b[611]), 99.0, 150, 0.99, "N_Protein", "Nucleoprotein", ref_seq_b, g_b["N_Protein"]),
    ]
    s3_cons_a = generate_consensus_sequence(ref_seq_a, s3_vars, "RSV-A_strain_Wuhan_2026", dropout=True, dropout_regions=[(15000, ref_length)])
    s3_cons_b = generate_consensus_sequence(ref_seq_b, s3_vars, "RSV-B_strain_Melbourne_2026", dropout=True, dropout_regions=[(14800, ref_length)])
    s3_cov_a_entry = _coverage_entry(212.4, 98.4, 97.2, 98.4, 63700, generate_amplicon_profile(ref_length, mean_depth=210), [[15000, ref_length]], ref_length, qc_thresholds)
    s3_cov_b_entry = _coverage_entry(164.2, 96.1, 94.8, 96.1, 49200, generate_amplicon_profile(ref_length, mean_depth=160), [[14800, ref_length]], ref_length, qc_thresholds)
    samples.append({
        "sample_id": "RSV-S003_COINFECTION",
        "qc_status": s3_cov_a_entry["qc_status"],
        "qc_reasons": s3_cov_a_entry["qc_reasons"],
        "primary_reference": "RSV-A_strain_Wuhan_2026",
        "qc": {
            "raw_reads": 154201,
            "raw_bases": 46260300,
            "trimmed_reads": 139502,
            "trimmed_bases": 41850600,
            "q30_rate_raw": 0.951,
            "q30_rate_trimmed": 0.972,
            "duplication_rate": 0.0512,
            "mean_read_length_raw": 150.0,
            "mean_read_length_trimmed": 149.1,
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": s3_cov_a_entry,
            "RSV-B_strain_Melbourne_2026": s3_cov_b_entry,
        },
        "variants": s3_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": _consensus_entry(ref_length, 48.2, s3_cons_a, [[15000, ref_length]]),
            "RSV-B_strain_Melbourne_2026": _consensus_entry(ref_length, 44.8, s3_cons_b, [[14800, ref_length]]),
        },
        "metadata": {"host": "Human", "collection_date": "2026-05-22", "location": "Brisbane"},
        "amplicons": [],
    })

    # Sample 4: Failed Run (FAIL on both)
    s4_cons_a = generate_consensus_sequence(ref_seq_a, [], "RSV-A_strain_Wuhan_2026", dropout=True, dropout_regions=[(1, ref_length)])
    s4_cons_b = generate_consensus_sequence(ref_seq_b, [], "RSV-B_strain_Melbourne_2026", dropout=True, dropout_regions=[(1, ref_length)])
    s4_cov_a_entry = _coverage_entry(4.1, 8.4, 1.1, 8.4, 1200, generate_amplicon_profile(ref_length, mean_depth=5, dropout=True, dropout_regions=[(1, ref_length)]), [[1, ref_length]], ref_length, qc_thresholds)
    s4_cov_b_entry = _coverage_entry(1.4, 2.1, 0.0, 2.1, 420, generate_amplicon_profile(ref_length, mean_depth=2, dropout=True, dropout_regions=[(1, ref_length)]), [[1, ref_length]], ref_length, qc_thresholds)
    samples.append({
        "sample_id": "RSV-S004_FAILED_RUN",
        "qc_status": s4_cov_a_entry["qc_status"],
        "qc_reasons": s4_cov_a_entry["qc_reasons"],
        "primary_reference": "RSV-A_strain_Wuhan_2026",
        "qc": {
            "raw_reads": 1402,
            "raw_bases": 420600,
            "trimmed_reads": 820,
            "trimmed_bases": 246000,
            "q30_rate_raw": 0.682,
            "q30_rate_trimmed": 0.714,
            "duplication_rate": 0.421,
            "mean_read_length_raw": 150.0,
            "mean_read_length_trimmed": 142.0,
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": s4_cov_a_entry,
            "RSV-B_strain_Melbourne_2026": s4_cov_b_entry,
        },
        "variants": [],
        "consensus": {
            "RSV-A_strain_Wuhan_2026": _consensus_entry(ref_length, 0.0, s4_cons_a, [[1, ref_length]]),
            "RSV-B_strain_Melbourne_2026": _consensus_entry(ref_length, 0.0, s4_cons_b, [[1, ref_length]]),
        },
        "metadata": {"host": "Control", "collection_date": "2026-05-22", "location": "Brisbane"},
        "amplicons": [],
    })

    # Compute mock pairwise distance matrices using consensus bases (null when zero mutual sites)
    distance_matrices = {}
    snp_distance_matrices = {}
    compared_sites_matrices = {}
    distance_matrices_by_name = {}
    sample_ids = [s["sample_id"] for s in samples]
    num_samples = len(samples)

    for ref in references:
        ref_name = ref["name"]
        ref_len = ref["length"]
        matrix = [[0.0] * num_samples for _ in range(num_samples)]
        snp_matrix = [[0] * num_samples for _ in range(num_samples)]
        sites_matrix = [[0] * num_samples for _ in range(num_samples)]

        for i in range(num_samples):
            seqA = samples[i]["consensus"][ref_name]["sequence"]
            sites_matrix[i][i] = sum(1 for c in seqA if c != "N")
            for j in range(i + 1, num_samples):
                seqB = samples[j]["consensus"][ref_name]["sequence"]
                diffs = 0
                compared = 0
                for pos_idx in range(ref_len):
                    charA = seqA[pos_idx]
                    charB = seqB[pos_idx]
                    if charA != "N" and charB != "N":
                        compared += 1
                        if charA.upper() != charB.upper():
                            diffs += 1
                dist = (diffs / compared) if compared > 0 else None
                snp_val = diffs if compared > 0 else None
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
        for i, sA_id in enumerate(sample_ids):
            ref_by_name[sA_id] = {}
            for j, sB_id in enumerate(sample_ids):
                ref_by_name[sA_id][sB_id] = matrix[i][j]
        distance_matrices_by_name[ref_name] = ref_by_name

    report = {
        "schema_version": parse_metrics.SCHEMA_VERSION,
        "run_id": "RSV_RUN_2026-05-30",
        "run_config": {
            "run_id": "RSV_RUN_2026-05-30",
            "platform": "illumina",
            "min_coverage": 10,
            "min_freq": 0.2,
            "min_consensus_depth": 10,
            "min_consensus_freq": 0.5,
            "min_map_qual": 20,
            "min_base_qual": 20,
            "min_variant_qual": 20,
            "use_iupac": False,
            "reference": "mock_rsv_a_b.fasta",
            "qc_thresholds": qc_thresholds,
        },
        "software_versions": {
            "bcftools": "1.25",
            "bwa": "0.7.19-r1273",
            "fastp": "1.4.0",
            "python": sys.version.split()[0],
            "samtools": "1.25",
        },
        "min_coverage_threshold": 10,
        "references": references,
        "genes": genes,
        "distance_matrices": distance_matrices,
        "distance_matrix_order": sample_ids,
        "distance_matrices_by_name": distance_matrices_by_name,
        "snp_distance_matrices": snp_distance_matrices,
        "compared_sites_matrices": compared_sites_matrices,
        "samples": samples,
    }

    parse_metrics.validate_report_schema(report)

    out_path = os.path.join(REPO_ROOT, "dashboard", "public", "mock_run_report.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"Multi-strain mock data successfully generated and validated:\n  - {out_path}")


if __name__ == "__main__":
    main()
