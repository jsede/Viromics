import json
import glob
import os
import sys

def resolve_results_dir():
    if len(sys.argv) > 1:
        return sys.argv[1]
    for candidate in ("test_ont_rsv_results", os.path.join("example_data", "ont_RSV", "results")):
        if os.path.exists(os.path.join(candidate, "viromics_run_report.json")):
            return candidate
    return "test_ont_rsv_results"

def main():
    results_dir = resolve_results_dir()
    report_path = os.path.join(results_dir, "viromics_run_report.json")
    
    # 1. Load run report JSON
    if not os.path.exists(report_path):
        print(f"Error: Report file {report_path} does not exist.", file=sys.stderr)
        sys.exit(1)
        
    try:
        with open(report_path) as f:
            r = json.load(f)
        print("Point 4 Check: Run report is valid JSON.")
    except Exception as e:
        print(f"Point 4 Fail: Failed to parse JSON: {e}", file=sys.stderr)
        sys.exit(1)
        
    # Check expected barcodes
    expected_barcodes = {"Sample13", "Sample29", "Sample35", "Sample40", "Sample45"}
    actual_barcodes = {s["sample_id"] for s in r["samples"]}
    if expected_barcodes.issubset(actual_barcodes):
        print(f"Point 4 Check: All expected barcodes {expected_barcodes} are present.")
    else:
        print(f"Point 4 Fail: Expected {expected_barcodes}, but got {actual_barcodes}", file=sys.stderr)
        sys.exit(1)
        
    # 2. Check coverage breadth relationship
    min_cov_threshold = r["min_coverage_threshold"]
    print(f"Min coverage threshold used: {min_cov_threshold}")
    
    point1_passed = True
    point2_passed = True
    
    reference_names = [ref["name"] for ref in r["references"]]
    print(f"References in report: {reference_names}")
    
    for s in r["samples"]:
        s_id = s["sample_id"]
        coverage = s["coverage"]
        consensus = s["consensus"]
        
        # Point 2 Check: Verify all references are explicitly present in coverage and consensus
        for ref_name in reference_names:
            if ref_name not in coverage:
                print(f"Point 2 Fail: Contig {ref_name} missing from coverage for sample {s_id}", file=sys.stderr)
                point2_passed = False
            if ref_name not in consensus:
                print(f"Point 2 Fail: Contig {ref_name} missing from consensus for sample {s_id}", file=sys.stderr)
                point2_passed = False
                
            if ref_name in coverage:
                cov_stats = coverage[ref_name]
                cov_10x = cov_stats["coverage_breadth_10x"]
                cov_min = cov_stats["coverage_breadth_min_cov"]
                
                # Point 1 Check: Verify cov_10x >= cov_min (since threshold 100 >= 10)
                if cov_10x < cov_min:
                    print(f"Point 1 Fail: Sample {s_id} Contig {ref_name} has cov_10x ({cov_10x}) < cov_min ({cov_min})", file=sys.stderr)
                    point1_passed = False
                    
                # Explicit representation of 0 depth if mapped_reads is 0
                if cov_stats.get("mapped_reads") == 0:
                    profile = cov_stats.get("coverage_profile", [])
                    all_zeros = all(p["depth"] == 0 for p in profile)
                    if not all_zeros:
                        print(f"Point 2 Fail: Sample {s_id} Contig {ref_name} has 0 mapped reads but non-zero coverage profile.", file=sys.stderr)
                        point2_passed = False
                        
    if point1_passed:
        print("Point 1 Check: coverage_breadth_10x >= coverage_breadth_min_cov holds true for all records.")
    if point2_passed:
        print("Point 2 Check: All reference contigs are explicitly represented for all samples (no omissions).")
        
    # 3. Check VCF genotypes
    vcf_paths = glob.glob(os.path.join(results_dir, "*", "variants", "*.vcf"))
    if not vcf_paths:
        print("Error: No VCF files found to validate genotypes.", file=sys.stderr)
        sys.exit(1)
        
    point3_passed = True
    for vcf in vcf_paths:
        with open(vcf) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                format_field = parts[8].split(":")
                sample_field = parts[9].split(":")
                
                if "GT" in format_field:
                    gt_idx = format_field.index("GT")
                    gt_val = sample_field[gt_idx]
                    if "/" in gt_val or "|" in gt_val:
                        print(f"Point 3 Fail: Found non-haploid genotype '{gt_val}' in {vcf} at line: {line.strip()}", file=sys.stderr)
                        point3_passed = False
                        
    if point3_passed:
        print("Point 3 Check: All VCF genotypes are strictly haploid (no '/' or '|').")
        
    if point1_passed and point2_passed and point3_passed:
        print("SUCCESS: All sanity checks passed perfectly!")
        sys.exit(0)
    else:
        print("FAILURE: One or more sanity checks failed.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
