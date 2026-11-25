#!/usr/bin/env python3
import glob
import json
import os
import sys

BIN_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(BIN_DIR, "..", ".."))
sys.path.insert(0, BIN_DIR)

import parse_metrics  # noqa: E402


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

    is_v2 = bool(r.get("schema_version"))
    if is_v2:
        parse_metrics.validate_report_schema(r)
        schema_path = os.path.join(REPO_ROOT, "schema", "viromics_run_report.schema.json")
        if os.path.exists(schema_path):
            try:
                import jsonschema
                with open(schema_path) as sf:
                    jsonschema.validate(instance=r, schema=json.load(sf))
                print("Schema v2 Check: Validated against viromics_run_report.schema.json.")
            except ImportError:
                print("Schema v2 Check: Validated via internal validate_report_schema (jsonschema pkg absent).")
    else:
        print("Notice: Legacy v1 report snapshot detected (missing schema_version); running structural sanity checks.")

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
                cov_min = cov_stats.get("coverage_breadth_min_cov", cov_stats.get(f"coverage_breadth_{min_cov_threshold}x"))
                if cov_min is None:
                    raise KeyError(f"Missing coverage_breadth_min_cov / coverage_breadth_{min_cov_threshold}x for {s_id}/{ref_name}")

                if min_cov_threshold >= 10 and cov_10x < cov_min:
                    print(
                        f"Point 1 Fail: Sample {s_id} Contig {ref_name} has cov_10x ({cov_10x}) < cov_min ({cov_min})",
                        file=sys.stderr,
                    )
                    point1_passed = False

                if cov_stats.get("mapped_reads") == 0:
                    profile = cov_stats.get("coverage_profile", [])
                    all_zeros = all(p["depth"] == 0 for p in profile)
                    if not all_zeros:
                        print(
                            f"Point 2 Fail: Sample {s_id} Contig {ref_name} has 0 mapped reads but non-zero coverage profile.",
                            file=sys.stderr,
                        )
                        point2_passed = False

    if point1_passed:
        print("Point 1 Check: coverage_breadth_10x >= coverage_breadth_min_cov holds true for all records.")
    if point2_passed:
        print("Point 2 Check: All reference contigs are explicitly represented for all samples (no omissions).")

    # 3. Check consensus/reported VCF genotypes
    consensus_vcfs = sorted(glob.glob(os.path.join(results_dir, "*", "variants", "*_consensus.vcf")))
    reported_vcfs = sorted(
        p for p in glob.glob(os.path.join(results_dir, "*", "variants", "*.vcf")) if not p.endswith("_consensus.vcf")
    )
    if not reported_vcfs:
        print("Error: No VCF files found to validate genotypes.", file=sys.stderr)
        sys.exit(1)

    use_iupac = bool(r.get("run_config", {}).get("use_iupac", False))
    point3_passed = True

    if consensus_vcfs:
        for vcf in consensus_vcfs:
            with open(vcf) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        continue
                    fmt = parts[8].split(":")
                    smp = parts[9].split(":")
                    if "GT" in fmt:
                        gt_val = smp[fmt.index("GT")]
                        if not use_iupac and ("/" in gt_val or "|" in gt_val):
                            print(
                                f"Point 3 Fail: Non-haploid GT '{gt_val}' in consensus VCF {vcf}: {line.strip()}",
                                file=sys.stderr,
                            )
                            point3_passed = False
        for vcf in reported_vcfs:
            with open(vcf) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    info = line.strip().split("\t")[7]
                    if "VAF=" not in info:
                        print(f"Point 3 Fail: Reported VCF {vcf} missing VAF tag: {line.strip()}", file=sys.stderr)
                        point3_passed = False
    else:
        for vcf in reported_vcfs:
            with open(vcf) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        continue
                    fmt = parts[8].split(":")
                    smp = parts[9].split(":")
                    if "GT" in fmt:
                        gt_val = smp[fmt.index("GT")]
                        if "/" in gt_val or "|" in gt_val:
                            print(f"Point 3 Fail: Non-haploid GT '{gt_val}' in {vcf}: {line.strip()}", file=sys.stderr)
                            point3_passed = False

    if point3_passed:
        print("Point 3 Check: VCF genotypes & decoupling invariants passed.")

    if point1_passed and point2_passed and point3_passed:
        print("SUCCESS: All sanity checks passed perfectly!")
        sys.exit(0)
    else:
        print("FAILURE: One or more sanity checks failed.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
