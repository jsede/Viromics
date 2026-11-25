#!/usr/bin/env python3
"""Regression tests for pairwise SNP distance matrices, deletion masking, compared_sites, and null zero-overlap."""

import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "bin"))
from parse_metrics import calculate_distance_matrices  # noqa: E402


def create_mock_data():
    ref_seq = "ACGT" * 1000  # 4000 bp
    ref_info = [{"name": "RefGenome", "length": len(ref_seq), "sequence": ref_seq}]

    samples_data = {
        "SampleA": {
            "sample_id": "SampleA",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": [],
        },
        "SampleB": {
            "sample_id": "SampleB",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": [
                {"reference": "RefGenome", "position": 10, "ref": "A", "alt": "T", "frequency": 1.0}
            ],
        },
        "SampleC": {
            "sample_id": "SampleC",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": [
                {"reference": "RefGenome", "position": 10, "ref": "A", "alt": "T", "frequency": 1.0},
                {"reference": "RefGenome", "position": 20, "ref": "C", "alt": "G", "frequency": 1.0},
            ],
        },
    }
    return samples_data, ref_info


def run_regression_test():
    print("Running stable distance matrix indexing regression test...")
    samples_data, ref_info = create_mock_data()

    matrices, order, matrices_by_name, snp_matrices, sites_matrices = calculate_distance_matrices(
        samples_data, ref_info, min_cov=10
    )

    assert order == ["SampleA", "SampleB", "SampleC"], f"Unexpected order: {order}"

    dist_ab = 1 / 4000
    dist_bc = 1 / 4000
    dist_ac = 2 / 4000

    matrix = matrices["RefGenome"]
    snps = snp_matrices["RefGenome"]
    sites = sites_matrices["RefGenome"]

    assert abs(matrix[0][1] - dist_ab) < 1e-9
    assert abs(matrix[1][2] - dist_bc) < 1e-9
    assert abs(matrix[0][2] - dist_ac) < 1e-9
    assert snps[0][1] == 1 and snps[1][2] == 1 and snps[0][2] == 2
    assert sites[0][0] == 4000 and sites[0][1] == 4000 and sites[0][2] == 4000

    by_name = matrices_by_name["RefGenome"]
    assert abs(by_name["SampleA"]["SampleB"] - dist_ab) < 1e-9
    assert abs(by_name["SampleB"]["SampleC"] - dist_bc) < 1e-9
    assert abs(by_name["SampleA"]["SampleC"] - dist_ac) < 1e-9

    # 2. Reordered samples preserve values when indexed by order
    reordered = {
        "SampleC": samples_data["SampleC"],
        "SampleA": samples_data["SampleA"],
        "SampleB": samples_data["SampleB"],
    }
    matrices_r, order_r, by_name_r, snp_r, sites_r = calculate_distance_matrices(
        reordered, ref_info, min_cov=10
    )
    assert order_r == ["SampleC", "SampleA", "SampleB"]
    idx_c, idx_a, idx_b = order_r.index("SampleC"), order_r.index("SampleA"), order_r.index("SampleB")
    assert abs(matrices_r["RefGenome"][idx_a][idx_b] - dist_ab) < 1e-9
    assert abs(matrices_r["RefGenome"][idx_b][idx_c] - dist_bc) < 1e-9
    assert abs(matrices_r["RefGenome"][idx_a][idx_c] - dist_ac) < 1e-9
    assert snp_r["RefGenome"][idx_a][idx_c] == 2
    assert sites_r["RefGenome"][idx_a][idx_c] == 4000
    assert abs(by_name_r["RefGenome"]["SampleA"]["SampleB"] - dist_ab) < 1e-9

    # 3. Deletion masking excludes deleted positions from both numerator and compared_sites
    samples_with_del = {
        "SampleA": samples_data["SampleA"],
        "SampleDel": {
            "sample_id": "SampleDel",
            "full_depths": {"RefGenome": [100] * 4000},
            "variants": [
                {"reference": "RefGenome", "position": 10, "ref": "A", "alt": "T", "frequency": 1.0},
                {"reference": "RefGenome", "position": 100, "ref": "GTCTGGTTTT", "alt": "G", "frequency": 1.0},
            ],
        },
    }
    _, _, by_name_del, snp_del, sites_del = calculate_distance_matrices(
        samples_with_del, ref_info, min_cov=10
    )
    assert sites_del["RefGenome"][0][1] == 3991
    assert snp_del["RefGenome"][0][1] == 1
    assert abs(by_name_del["RefGenome"]["SampleA"]["SampleDel"] - (1 / 3991)) < 1e-9

    # 4. Zero-overlap yields None (JSON null) for distance & SNP count and 0 for compared_sites
    samples_zero = {
        "SampleA": samples_data["SampleA"],
        "SampleUncovered": {
            "sample_id": "SampleUncovered",
            "full_depths": {"RefGenome": [0] * 4000},
            "variants": [],
        },
    }
    _, _, by_name_zero, snp_zero, sites_zero = calculate_distance_matrices(
        samples_zero, ref_info, min_cov=10
    )
    assert by_name_zero["RefGenome"]["SampleA"]["SampleUncovered"] is None
    assert snp_zero["RefGenome"][0][1] is None
    assert sites_zero["RefGenome"][0][1] == 0
    assert by_name_zero["RefGenome"]["SampleUncovered"]["SampleUncovered"] == 0.0
    assert snp_zero["RefGenome"][1][1] == 0

    print("SUCCESS: Regression test passed (deletion masking, compared_sites, and null zero-overlap verified).")


if __name__ == "__main__":
    run_regression_test()
