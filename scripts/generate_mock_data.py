#!/usr/bin/env python3
import os
import json
import random
import math

def generate_amplicon_profile(length, num_bins=500, mean_depth=300, dropout=False, dropout_regions=None):
    """Simulates realistic amplicon coverage depth with peaks and troughs using sine and random noise."""
    profile = []
    bin_size = length / num_bins
    
    num_amplicons = 20
    amp_positions = [i * (length / num_amplicons) for i in range(num_amplicons)]
    amp_heights = [random.randint(int(mean_depth*0.5), int(mean_depth*1.5)) for _ in range(num_amplicons)]
    
    if dropout and dropout_regions:
        for idx in range(num_amplicons):
            pos = amp_positions[idx]
            for start, end in dropout_regions:
                if start <= pos <= end:
                    amp_heights[idx] = random.randint(0, 3)
    
    for i in range(num_bins):
        pos = i * bin_size
        depth = 2 # baseline noise
        for amp_pos, amp_height in zip(amp_positions, amp_heights):
            dist = pos - amp_pos
            width = length / 30
            weight = math.exp(-(dist**2) / (2 * (width**2)))
            depth += amp_height * weight
            
        depth = max(1, depth * random.uniform(0.85, 1.15))
        if dropout and dropout_regions:
            for start, end in dropout_regions:
                if start <= pos <= end:
                    depth = max(0, depth * random.uniform(0.01, 0.05))
                    
        profile.append({
            "pos": int(pos) + 1,
            "depth": round(depth)
        })
    return profile

def generate_consensus_sequence(ref_seq, variants, ref_name, dropout=False, dropout_regions=None):
    seq_list = list(ref_seq)
    for var in variants:
        if var["reference"] == ref_name:
            pos_idx = var["position"] - 1
            if 0 <= pos_idx < len(seq_list):
                seq_list[pos_idx] = var["alt"]
    if dropout and dropout_regions:
        for start, end in dropout_regions:
            for i in range(start - 1, end):
                if 0 <= i < len(seq_list):
                    seq_list[i] = 'N'
    return "".join(seq_list)

def main():
    print("Generating Multi-Strain RSV-A/B Mock Run Report...")
    random.seed(42)
    
    # 1. References list
    ref_length = 15200
    bases = ['A', 'C', 'G', 'T']
    ref_seq_a = "".join(random.choices(bases, k=ref_length))
    ref_seq_b = "".join(random.choices(bases, k=ref_length))
    
    references = [
        {"name": "RSV-A_strain_Wuhan_2026", "length": ref_length, "sequence": ref_seq_a},
        {"name": "RSV-B_strain_Melbourne_2026", "length": ref_length, "sequence": ref_seq_b}
    ]
    
    # 2. Genes grouped by reference
    genes = {
        "RSV-A_strain_Wuhan_2026": [
            {"name": "N_Protein", "type": "gene", "start": 100, "end": 1200, "strand": "+", "product": "Nucleoprotein"},
            {"name": "F_Protein", "type": "gene", "start": 1500, "end": 3500, "strand": "+", "product": "Fusion glycoprotein"},
            {"name": "G_Protein", "type": "gene", "start": 4000, "end": 6000, "strand": "+", "product": "Attachment glycoprotein"},
            {"name": "L_Protein", "type": "gene", "start": 7000, "end": 14500, "strand": "+", "product": "RNA-directed RNA polymerase"}
        ],
        "RSV-B_strain_Melbourne_2026": [
            {"name": "N_Protein", "type": "gene", "start": 100, "end": 1200, "strand": "+", "product": "Nucleoprotein"},
            {"name": "F_Protein", "type": "gene", "start": 1500, "end": 3500, "strand": "+", "product": "Fusion glycoprotein"},
            {"name": "G_Protein", "type": "gene", "start": 3800, "end": 5800, "strand": "+", "product": "Attachment glycoprotein"},
            {"name": "L_Protein", "type": "gene", "start": 7000, "end": 14500, "strand": "+", "product": "RNA-directed RNA polymerase"}
        ]
    }
    
    # 3. Create Sample Profiles
    samples = []
    
    # Sample 1: Pure RSV-A Infection (PASS on RSV-A, FAIL on RSV-B)
    s1_cov_a = generate_amplicon_profile(ref_length, mean_depth=350)
    s1_cov_b = generate_amplicon_profile(ref_length, mean_depth=1)
    
    s1_vars = [
        {"reference": "RSV-A_strain_Wuhan_2026", "position": 520, "ref": "C", "alt": "T", "quality": 99.0, "depth": 310, "frequency": 0.99, "gene": "N_Protein", "type": "Coding"},
        {"reference": "RSV-A_strain_Wuhan_2026", "position": 2105, "ref": "A", "alt": "G", "quality": 95.0, "depth": 280, "frequency": 0.97, "gene": "F_Protein", "type": "Coding"},
        {"reference": "RSV-A_strain_Wuhan_2026", "position": 4820, "ref": "G", "alt": "C", "quality": 88.0, "depth": 410, "frequency": 0.94, "gene": "G_Protein", "type": "Coding"}
    ]
    
    s1_cons_a = generate_consensus_sequence(ref_seq_a, s1_vars, "RSV-A_strain_Wuhan_2026")
    s1_cons_b = generate_consensus_sequence(ref_seq_b, s1_vars, "RSV-B_strain_Melbourne_2026", dropout=True, dropout_regions=[(1, ref_length)])
    
    samples.append({
        "sample_id": "RSV-S001_A_Positive",
        "qc": {
            "raw_reads": 128042,
            "raw_bases": 38412600,
            "trimmed_reads": 112040,
            "trimmed_bases": 33612000,
            "q30_rate_raw": 94.8,
            "q30_rate_trimmed": 96.9,
            "duplication_rate": 6.42
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": {
                "mean_depth": 352.4,
                "coverage_breadth_10x": 99.80,
                "coverage_breadth_30x": 99.20,
                "coverage_profile": s1_cov_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "mean_depth": 0.8,
                "coverage_breadth_10x": 0.00,
                "coverage_breadth_30x": 0.00,
                "coverage_profile": s1_cov_b
            }
        },
        "variants": s1_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": {
                "length": ref_length,
                "gc_content": 48.2,
                "n_content": 0.20,
                "sequence": s1_cons_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "length": ref_length,
                "gc_content": 45.1,
                "n_content": 100.00,
                "sequence": s1_cons_b
            }
        }
    })
    
    # Sample 2: Pure RSV-B Infection (PASS on RSV-B, FAIL on RSV-A)
    s2_cov_a = generate_amplicon_profile(ref_length, mean_depth=1)
    s2_cov_b = generate_amplicon_profile(ref_length, mean_depth=280)
    
    s2_vars = [
        {"reference": "RSV-B_strain_Melbourne_2026", "position": 612, "ref": "T", "alt": "C", "quality": 99.0, "depth": 240, "frequency": 0.99, "gene": "N_Protein", "type": "Coding"},
        {"reference": "RSV-B_strain_Melbourne_2026", "position": 4105, "ref": "G", "alt": "A", "quality": 92.0, "depth": 320, "frequency": 0.95, "gene": "G_Protein", "type": "Coding"}
    ]
    
    s2_cons_a = generate_consensus_sequence(ref_seq_a, s2_vars, "RSV-A_strain_Wuhan_2026", dropout=True, dropout_regions=[(1, ref_length)])
    s2_cons_b = generate_consensus_sequence(ref_seq_b, s2_vars, "RSV-B_strain_Melbourne_2026")
    
    samples.append({
        "sample_id": "RSV-S002_B_Positive",
        "qc": {
            "raw_reads": 98402,
            "raw_bases": 29520600,
            "trimmed_reads": 86240,
            "trimmed_bases": 25872000,
            "q30_rate_raw": 93.9,
            "q30_rate_trimmed": 96.1,
            "duplication_rate": 8.21
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": {
                "mean_depth": 0.6,
                "coverage_breadth_10x": 0.00,
                "coverage_breadth_30x": 0.00,
                "coverage_profile": s2_cov_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "mean_depth": 281.2,
                "coverage_breadth_10x": 99.50,
                "coverage_breadth_30x": 98.40,
                "coverage_profile": s2_cov_b
            }
        },
        "variants": s2_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": {
                "length": ref_length,
                "gc_content": 48.2,
                "n_content": 100.00,
                "sequence": s2_cons_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "length": ref_length,
                "gc_content": 44.8,
                "n_content": 0.50,
                "sequence": s2_cons_b
            }
        }
    })
    
    # Sample 3: RSV-A and RSV-B Co-infection / Co-detection! (PASS on both!)
    s3_cov_a = generate_amplicon_profile(ref_length, mean_depth=210)
    s3_cov_b = generate_amplicon_profile(ref_length, mean_depth=160)
    
    s3_vars = [
        {"reference": "RSV-A_strain_Wuhan_2026", "position": 520, "ref": "C", "alt": "T", "quality": 99.0, "depth": 180, "frequency": 0.99, "gene": "N_Protein", "type": "Coding"},
        {"reference": "RSV-B_strain_Melbourne_2026", "position": 612, "ref": "T", "alt": "C", "quality": 99.0, "depth": 150, "frequency": 0.99, "gene": "N_Protein", "type": "Coding"}
    ]
    
    s3_cons_a = generate_consensus_sequence(ref_seq_a, s3_vars, "RSV-A_strain_Wuhan_2026")
    s3_cons_b = generate_consensus_sequence(ref_seq_b, s3_vars, "RSV-B_strain_Melbourne_2026")
    
    samples.append({
        "sample_id": "RSV-S003_COINFECTION",
        "qc": {
            "raw_reads": 154201,
            "raw_bases": 46260300,
            "trimmed_reads": 139502,
            "trimmed_bases": 41850600,
            "q30_rate_raw": 95.1,
            "q30_rate_trimmed": 97.2,
            "duplication_rate": 5.12
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": {
                "mean_depth": 212.4,
                "coverage_breadth_10x": 98.40,
                "coverage_breadth_30x": 97.20,
                "coverage_profile": s3_cov_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "mean_depth": 164.2,
                "coverage_breadth_10x": 96.10,
                "coverage_breadth_30x": 94.80,
                "coverage_profile": s3_cov_b
            }
        },
        "variants": s3_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": {
                "length": ref_length,
                "gc_content": 48.2,
                "n_content": 1.60,
                "sequence": s3_cons_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "length": ref_length,
                "gc_content": 44.8,
                "n_content": 3.90,
                "sequence": s3_cons_b
            }
        }
    })
    
    # Sample 4: Low Quality Fail (reads < 2000, coverage < 10%)
    s4_cov_a = generate_amplicon_profile(ref_length, mean_depth=5, dropout=True, dropout_regions=[(1, ref_length)])
    s4_cov_b = generate_amplicon_profile(ref_length, mean_depth=2, dropout=True, dropout_regions=[(1, ref_length)])
    s4_vars = []
    s4_cons_a = generate_consensus_sequence(ref_seq_a, s4_vars, "RSV-A_strain_Wuhan_2026", dropout=True, dropout_regions=[(1, ref_length)])
    s4_cons_b = generate_consensus_sequence(ref_seq_b, s4_vars, "RSV-B_strain_Melbourne_2026", dropout=True, dropout_regions=[(1, ref_length)])
    
    samples.append({
        "sample_id": "RSV-S004_FAILED_RUN",
        "qc": {
            "raw_reads": 1402,
            "raw_bases": 420600,
            "trimmed_reads": 820,
            "trimmed_bases": 246000,
            "q30_rate_raw": 68.2,
            "q30_rate_trimmed": 71.4,
            "duplication_rate": 42.10
        },
        "coverage": {
            "RSV-A_strain_Wuhan_2026": {
                "mean_depth": 4.1,
                "coverage_breadth_10x": 8.40,
                "coverage_breadth_30x": 1.10,
                "coverage_profile": s4_cov_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "mean_depth": 1.4,
                "coverage_breadth_10x": 2.10,
                "coverage_breadth_30x": 0.00,
                "coverage_profile": s4_cov_b
            }
        },
        "variants": s4_vars,
        "consensus": {
            "RSV-A_strain_Wuhan_2026": {
                "length": ref_length,
                "gc_content": 48.2,
                "n_content": 100.00,
                "sequence": s4_cons_a
            },
            "RSV-B_strain_Melbourne_2026": {
                "length": ref_length,
                "gc_content": 44.8,
                "n_content": 100.00,
                "sequence": s4_cons_b
            }
        }
    })
    
    # 4. Save report
    report = {
        "run_id": "RSV_RUN_2026-05-30",
        "references": references,
        "genes": genes,
        "samples": samples
    }
    
    os.makedirs("dashboard/public", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    
    with open("dashboard/public/mock_run_report.json", 'w') as f:
        json.dump(report, f, indent=2)
        
    with open("results/mock_run_report.json", 'w') as f:
        json.dump(report, f, indent=2)
        
    print("Multi-strain mock data successfully generated:")
    print("  - dashboard/public/mock_run_report.json")
    print("  - results/mock_run_report.json")

if __name__ == "__main__":
    main()
