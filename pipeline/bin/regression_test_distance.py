#!/usr/bin/env python3
import sys
import os
import json

# Add pipeline/bin to path to import calculate_distance_matrices
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "bin"))
try:
    from parse_metrics import calculate_distance_matrices
except ImportError:
    # Fallback to local import if run from within pipeline/bin
    from parse_metrics import calculate_distance_matrices

def create_mock_data():
    # Define a reference genome
    ref_seq = "ACGT" * 1000  # 4000 bp
    ref_info = [{"name": "RefGenome", "length": len(ref_seq), "sequence": ref_seq}]
    
    # Define 3 samples with distinct SNP profiles
    # Sample A: no SNPs
    # Sample B: SNP at pos 10 (A->T)
    # Sample C: SNP at pos 10 (A->T) and pos 20 (C->G)
    
    samples_data = {
        "SampleA": {
            "sample_id": "SampleA",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": []
        },
        "SampleB": {
            "sample_id": "SampleB",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": [
                {"reference": "RefGenome", "position": 10, "ref": "A", "alt": "T"}
            ]
        },
        "SampleC": {
            "sample_id": "SampleC",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": [
                {"reference": "RefGenome", "position": 10, "ref": "A", "alt": "T"},
                {"reference": "RefGenome", "position": 20, "ref": "C", "alt": "G"}
            ]
        }
    }
    return samples_data, ref_info

def run_regression_test():
    print("Running stable distance matrix indexing regression test...")
    samples_data, ref_info = create_mock_data()
    ref_seq = ref_info[0]["sequence"]
    
    # 1. Test standard order
    matrices, order, matrices_by_name = calculate_distance_matrices(samples_data, ref_info, min_cov=10)
    
    assert order == ["SampleA", "SampleB", "SampleC"], f"Expected order to be ['SampleA', 'SampleB', 'SampleC'], got {order}"
    
    # Verify distances manually
    # Sample A vs Sample B: 1 difference at position 10 (A vs T). 4000 sites compared.
    # Dist = 1 / 4000 = 0.00025
    dist_ab = 1 / 4000
    # Sample B vs Sample C: 1 difference at position 20 (C vs G). 4000 sites compared.
    # Dist = 1 / 4000 = 0.00025
    dist_bc = 1 / 4000
    # Sample A vs Sample C: 2 differences at pos 10 and 20.
    # Dist = 2 / 4000 = 0.0005
    dist_ac = 2 / 4000
    
    # Check index based matrix
    matrix = matrices["RefGenome"]
    assert abs(matrix[0][1] - dist_ab) < 1e-9, f"Expected A vs B dist to be {dist_ab}, got {matrix[0][1]}"
    assert abs(matrix[1][2] - dist_bc) < 1e-9, f"Expected B vs C dist to be {dist_bc}, got {matrix[1][2]}"
    assert abs(matrix[0][2] - dist_ac) < 1e-9, f"Expected A vs C dist to be {dist_ac}, got {matrix[0][2]}"
    
    # Check name based matrix
    by_name = matrices_by_name["RefGenome"]
    assert abs(by_name["SampleA"]["SampleB"] - dist_ab) < 1e-9
    assert abs(by_name["SampleB"]["SampleC"] - dist_bc) < 1e-9
    assert abs(by_name["SampleA"]["SampleC"] - dist_ac) < 1e-9
    
    # 2. Test reordered samples (to simulate Nextflow non-deterministic input sorting)
    reordered_samples_data = {
        "SampleC": samples_data["SampleC"],
        "SampleA": samples_data["SampleA"],
        "SampleB": samples_data["SampleB"]
    }
    
    matrices_r, order_r, matrices_by_name_r = calculate_distance_matrices(reordered_samples_data, ref_info, min_cov=10)
    
    assert order_r == ["SampleC", "SampleA", "SampleB"], f"Expected order to be ['SampleC', 'SampleA', 'SampleB'], got {order_r}"
    
    matrix_r = matrices_r["RefGenome"]
    
    # Using order_r to map sample names to their corresponding matrix indices
    idx_c = order_r.index("SampleC")
    idx_a = order_r.index("SampleA")
    idx_b = order_r.index("SampleB")
    
    assert idx_c == 0
    assert idx_a == 1
    assert idx_b == 2
    
    # The matrix indices should match the new order, but lookup using the names and index mapping should still yield the correct values
    assert abs(matrix_r[idx_a][idx_b] - dist_ab) < 1e-9, "Alignment failed for A vs B after reordering"
    assert abs(matrix_r[idx_b][idx_c] - dist_bc) < 1e-9, "Alignment failed for B vs C after reordering"
    assert abs(matrix_r[idx_a][idx_c] - dist_ac) < 1e-9, "Alignment failed for A vs C after reordering"
    
    # The name-based map should also remain completely identical regardless of insertion order
    by_name_r = matrices_by_name_r["RefGenome"]
    assert abs(by_name_r["SampleA"]["SampleB"] - dist_ab) < 1e-9
    assert abs(by_name_r["SampleB"]["SampleC"] - dist_bc) < 1e-9
    assert abs(by_name_r["SampleA"]["SampleC"] - dist_ac) < 1e-9
    
    # 3. Test deletion masking (deletions excluded from SNP distance comparisons)
    samples_with_del = {
        "SampleA": samples_data["SampleA"],
        "SampleDel": {
            "sample_id": "SampleDel",
            "full_depths": {"RefGenome": [100] * len(ref_seq)},
            "variants": [
                # SNP at pos 10
                {"reference": "RefGenome", "position": 10, "ref": "A", "alt": "T", "frequency": 1.0},
                # 9-bp deletion at pos 100 (bases 101..109 deleted)
                {"reference": "RefGenome", "position": 100, "ref": "GTCTGGTTTT", "alt": "G", "frequency": 1.0}
            ]
        }
    }
    matrices_del, _, by_name_del = calculate_distance_matrices(samples_with_del, ref_info, min_cov=10)
    # Total sites compared should exclude the 9 deleted positions (4000 - 9 = 3991)
    # Number of SNP differences = 1 (at pos 10)
    expected_dist_del = 1 / (4000 - 9)
    assert abs(by_name_del["RefGenome"]["SampleA"]["SampleDel"] - expected_dist_del) < 1e-9, f"Expected deletion-aware dist {expected_dist_del}, got {by_name_del['RefGenome']['SampleA']['SampleDel']}"

    # 4. Test zero-overlap handling (unsequenced sample has dist 1.0 to other samples, 0.0 to itself)
    samples_with_zero_cov = {
        "SampleA": samples_data["SampleA"],
        "SampleUncovered": {
            "sample_id": "SampleUncovered",
            "full_depths": {"RefGenome": [0] * len(ref_seq)},
            "variants": []
        }
    }
    _, _, by_name_zero = calculate_distance_matrices(samples_with_zero_cov, ref_info, min_cov=10)
    assert by_name_zero["RefGenome"]["SampleA"]["SampleUncovered"] == 1.0, f"Expected 1.0 for zero overlap, got {by_name_zero['RefGenome']['SampleA']['SampleUncovered']}"
    assert by_name_zero["RefGenome"]["SampleUncovered"]["SampleUncovered"] == 0.0, f"Expected 0.0 for self-distance, got {by_name_zero['RefGenome']['SampleUncovered']['SampleUncovered']}"

    print("SUCCESS: Regression test passed perfectly! Sample names align correctly with matrix rows/columns regardless of order.")

if __name__ == "__main__":
    run_regression_test()
