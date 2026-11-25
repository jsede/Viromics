#!/usr/bin/env python3
"""
End-to-end biological & schema verification script for Viromics CI synthetic fixture runs.

Validates:
  1. JSON schema conformance (both internal parse_metrics.validate_report_schema & jsonschema).
  2. Provenance & tool versions (run_config, software_versions).
  3. SynthPos: PASS QC tier, 6/6 passing amplicons, clonal SNP 200 G>A (ORF1: A34T),
     ~31% subclonal insertion 370 T>TA reported in Illumina VCF (is_consensus=False)
     and excluded from consensus in BOTH Illumina & ONT, clonal 6-bp in-frame deletion
     at 571 preserved cleanly (consensus len 1424 bp), and 0 variant at primer-decoy pos 10.
  4. BAM alignment audit via samtools view on SynthPos & SynthDropout trimmed BAMs:
     0% of primary mapped reads (-F 0x904) overlap any primer interval in primer BED.
  5. SynthDropout: internal dropout gap covering 716..945 in low_coverage_intervals and
     N-masked in consensus, Amp 4 status == 'fail' while Amps 1,2,3,5,6 == 'pass'.
  6. Distance matrices: SynthPos vs SynthDropout snp_distance == 1 with >=1000 compared sites;
     SynthNTC pairs have compared_sites == 0 and distance is None.
  7. SynthNTC: empty 0-read NTC completes cleanly, raw_reads == 0, qc_status == 'FAIL',
     100% consensus N-content.
"""

import argparse
import json
import os
import re
import subprocess
import sys

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO_ROOT, "pipeline", "bin"))

import parse_metrics  # noqa: E402


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
REF_CONSUMING_OPS = {"M", "D", "N", "=", "X"}


def load_primer_intervals(primer_bed):
    intervals = []
    with open(primer_bed, "r") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            strand = parts[5] if len(parts) > 5 else ("-" if "RIGHT" in parts[3].upper() else "+")
            intervals.append((parts[0], int(parts[1]), int(parts[2]), parts[3], strand))
    return intervals


def cigar_ref_len(cigar):
    return sum(int(n) for n, op in CIGAR_RE.findall(cigar) if op in REF_CONSUMING_OPS)


def check_bam_zero_primer_overlap(bam_path, primers, ref_name="SYN_VIRUS_1"):
    cmd = ["samtools", "view", "-F", "0x904", bam_path]
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    overlapping = []
    total_reads = 0
    left_primers = [p for p in primers if p[0] == ref_name and p[4] == "+"]
    right_primers = [p for p in primers if p[0] == ref_name and p[4] == "-"]
    min_insert_start = min(p[2] for p in left_primers) if left_primers else 0
    max_insert_end = max(p[1] for p in right_primers) if right_primers else 10**9

    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        total_reads += 1
        cols = line.split("\t")
        rname = cols[2]
        if rname != ref_name:
            continue
        pos0 = int(cols[3]) - 1
        rlen = cigar_ref_len(cols[5])
        read_start = pos0
        read_end = pos0 + rlen

        # 1. Exterior terminal primers must never have residual aligned bases
        if read_start < min_insert_start or read_end > max_insert_end:
            overlapping.append((cols[0], read_start, read_end, "exterior_primer", min_insert_start, max_insert_end))
            continue

        # 2. Internal tiled amplicons: aligned read start must not sit inside a '+' primer,
        # and aligned read end must not sit inside a '-' primer
        for _c, p_s, p_e, p_name, _st in left_primers:
            if p_s <= read_start < p_e:
                overlapping.append((cols[0], read_start, read_end, p_name, p_s, p_e))
                break
        else:
            for _c, p_s, p_e, p_name, _st in right_primers:
                if p_s < read_end <= p_e:
                    overlapping.append((cols[0], read_start, read_end, p_name, p_s, p_e))
                    break

    assert total_reads > 0, f"Expected mapped reads in {bam_path}, found 0."
    assert not overlapping, (
        f"Residual primer boundary overlap in {bam_path}: {len(overlapping)}/{total_reads} reads "
        f"remain inside primer footprints (first example: {overlapping[0]})"
    )
    return total_reads


def interval_covers_range(intervals, target_start, target_end):
    for iv_s, iv_e in intervals:
        if iv_s <= target_start and iv_e >= target_end:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description="Verify synthetic Viromics E2E run outputs.")
    parser.add_argument("--outdir", required=True, help="Pipeline output directory")
    parser.add_argument("--platform", required=True, choices=["illumina", "ont"], help="Platform run mode")
    parser.add_argument(
        "--primer-bed",
        default=os.path.join(REPO_ROOT, "pipeline", "test_data", "primers_synth.bed"),
        help="Path to primer BED",
    )
    parser.add_argument(
        "--schema",
        default=os.path.join(REPO_ROOT, "schema", "viromics_run_report.schema.json"),
        help="Path to report JSON schema",
    )
    parser.add_argument("--require-jsonschema", action="store_true", help="Fail if jsonschema package is absent")
    args = parser.parse_args()

    report_path = os.path.join(args.outdir, "viromics_run_report.json")
    assert os.path.exists(report_path), f"Missing run report: {report_path}"
    with open(report_path, "r") as f:
        report = json.load(f)
    with open(args.schema, "r") as f:
        schema = json.load(f)

    # 1. Schema validation
    parse_metrics.validate_report_schema(report)
    try:
        import jsonschema
        jsonschema.validate(instance=report, schema=schema)
    except ImportError:
        if args.require_jsonschema:
            raise RuntimeError("jsonschema package is required (--require-jsonschema) but not installed.")

    # 2. Provenance assertions
    assert report["run_config"]["platform"] == args.platform
    assert len(report["software_versions"]) >= 3, f"Unexpectedly sparse software_versions: {report['software_versions']}"
    ref_name = "SYN_VIRUS_1"
    samples = {s["sample_id"]: s for s in report["samples"]}
    assert set(samples.keys()) == {"SynthPos", "SynthDropout", "SynthNTC"}, f"Unexpected samples: {list(samples.keys())}"

    # 3. SynthPos assertions
    pos_s = samples["SynthPos"]
    assert pos_s["qc_status"] == "PASS", f"SynthPos expected PASS, got {pos_s['qc_status']} ({pos_s['qc_reasons']})"
    assert pos_s["primary_reference"] == ref_name
    pos_amps = {a["amplicon_id"]: a for a in pos_s["amplicons"]}
    assert len(pos_amps) == 6, f"Expected 6 amplicons for SynthPos, got {len(pos_amps)}"
    for aid, a in pos_amps.items():
        assert a["status"] == "pass", f"SynthPos amplicon {aid} expected pass, got {a}"

    pos_variants = {(v["position"], v["ref"], v["alt"]): v for v in pos_s["variants"]}
    # Decoy pos 10 inside Synth_1_LEFT MUST NOT appear in variants
    assert all(v["position"] > 25 for v in pos_s["variants"]), (
        f"Residual primer variant detected <=25 bp: {[v for v in pos_s['variants'] if v['position'] <= 25]}"
    )

    # Clonal SNP 200 G>A -> ORF1: A34T
    assert (200, "G", "A") in pos_variants, f"Missing clonal SNP 200 G>A in SynthPos variants: {list(pos_variants.keys())}"
    snp200 = pos_variants[(200, "G", "A")]
    assert snp200["is_consensus"] is True
    assert snp200["aa_change"] == "A34T", f"Expected A34T for pos 200, got {snp200}"
    assert snp200["consequences"][0]["gene"] == "ORF1"
    assert snp200["consequences"][0]["effect"] == "missense"

    # Clonal 6-bp in-frame deletion at anchor 571 (CATGCTG -> C)
    assert (571, "CATGCTG", "C") in pos_variants, f"Missing 6-bp deletion at 571 in SynthPos variants: {list(pos_variants.keys())}"
    del571 = pos_variants[(571, "CATGCTG", "C")]
    assert del571["is_consensus"] is True
    assert del571["consequences"][0]["effect"] == "inframe_deletion"

    # Subclonal insertion 370 T>TA (~31% AF)
    if args.platform == "illumina":
        assert (370, "T", "TA") in pos_variants, f"Missing subclonal 370 T>TA in Illumina variants: {list(pos_variants.keys())}"
        ins370 = pos_variants[(370, "T", "TA")]
        assert ins370["is_consensus"] is False, f"Subclonal 370 T>TA must NOT be consensus: {ins370}"
        assert 0.20 <= ins370["frequency"] < 0.50, f"Unexpected frequency for 370 T>TA: {ins370['frequency']}"

    # Consensus FASTA preserves clonal 6-bp del (1424 bp) & clonal SNP, excludes 370 +A insertion, zero N across deletion
    pos_cons = pos_s["consensus"][ref_name]
    assert pos_cons["length"] == 1424, f"SynthPos consensus expected 1424 bp (-6 bp del, no +1 ins), got {pos_cons['length']}"
    assert pos_cons["internal_n_count"] == 0, f"SynthPos expected 0 internal Ns, got {pos_cons['internal_n_count']}"

    # Verify consensus VCF excludes 370 T>TA
    cons_vcf = os.path.join(args.outdir, "SynthPos", "variants", "SynthPos_consensus.vcf")
    assert os.path.exists(cons_vcf), f"Missing {cons_vcf}"
    with open(cons_vcf, "r") as f:
        cons_vcf_positions = {int(l.split("\t")[1]) for l in f if l.strip() and not l.startswith("#")}
    assert 200 in cons_vcf_positions and 571 in cons_vcf_positions, f"Missing clonal sites in consensus VCF: {cons_vcf_positions}"
    assert 370 not in cons_vcf_positions, f"Subclonal indel 370 leaked into consensus VCF: {cons_vcf_positions}"

    # 4. Direct BAM primer-clipping verification (0% residual primer overlap)
    primers = load_primer_intervals(args.primer_bed)
    pos_reads = check_bam_zero_primer_overlap(
        os.path.join(args.outdir, "SynthPos", "alignment", "SynthPos_trimmed.bam"), primers, ref_name
    )
    drop_reads = check_bam_zero_primer_overlap(
        os.path.join(args.outdir, "SynthDropout", "alignment", "SynthDropout_trimmed.bam"), primers, ref_name
    )

    # 5. SynthDropout assertions (Amp 4 gap 716..945 masked as N & reported in low_coverage_intervals)
    drop_s = samples["SynthDropout"]
    drop_intervals = drop_s["coverage"][ref_name]["low_coverage_intervals"]
    assert interval_covers_range(drop_intervals, 716, 945), (
        f"SynthDropout low_coverage_intervals {drop_intervals} failed to cover dropout 716..945"
    )
    assert drop_s["consensus"][ref_name]["internal_n_count"] >= 220, (
        f"SynthDropout expected >=220 internal Ns across 230-bp gap, got {drop_s['consensus'][ref_name]['internal_n_count']}"
    )
    drop_amps = {a["amplicon_id"]: a for a in drop_s["amplicons"]}
    assert drop_amps["Synth_4"]["status"] == "fail", f"Expected Synth_4 to fail in SynthDropout: {drop_amps['Synth_4']}"
    for pass_amp in ("Synth_1", "Synth_2", "Synth_3", "Synth_5", "Synth_6"):
        assert drop_amps[pass_amp]["status"] == "pass", f"Expected {pass_amp} pass in SynthDropout: {drop_amps[pass_amp]}"

    # 6. SynthNTC assertions (0-read sample finishes gracefully and classifies FAIL)
    ntc_s = samples["SynthNTC"]
    assert ntc_s["qc"]["raw_reads"] == 0, f"SynthNTC expected 0 raw_reads, got {ntc_s['qc']}"
    assert ntc_s["qc_status"] == "FAIL", f"SynthNTC expected FAIL, got {ntc_s['qc_status']}"
    assert ntc_s["consensus"][ref_name]["n_content"] == 100.0

    # 7. Distance matrix assertions
    order = report["distance_matrix_order"]
    i_pos = order.index("SynthPos")
    i_drop = order.index("SynthDropout")
    i_ntc = order.index("SynthNTC")

    snp_mat = report["snp_distance_matrices"][ref_name]
    sites_mat = report["compared_sites_matrices"][ref_name]
    dist_mat = report["distance_matrices"][ref_name]

    assert snp_mat[i_pos][i_drop] == 1, f"Expected 1 SNP diff (pos 200) SynthPos vs SynthDropout, got {snp_mat[i_pos][i_drop]}"
    assert sites_mat[i_pos][i_drop] >= 1000, f"Expected >=1000 compared sites, got {sites_mat[i_pos][i_drop]}"
    assert dist_mat[i_pos][i_drop] is not None and dist_mat[i_pos][i_drop] > 0

    assert sites_mat[i_pos][i_ntc] == 0, f"Expected 0 compared sites vs NTC, got {sites_mat[i_pos][i_ntc]}"
    assert dist_mat[i_pos][i_ntc] is None, f"Expected None distance vs NTC, got {dist_mat[i_pos][i_ntc]}"
    assert snp_mat[i_pos][i_ntc] is None, f"Expected None SNP count vs NTC, got {snp_mat[i_pos][i_ntc]}"

    print(
        f"SUCCESS [{args.platform}]: Schema + all biological ground-truth assertions passed "
        f"(BAM primary reads verified 0% primer overlap: SynthPos={pos_reads}, SynthDropout={drop_reads})."
    )


if __name__ == "__main__":
    main()
