#!/usr/bin/env python3
"""
Decouples reported subclonal/major variants from consensus variant selection.

- Reported VCF keeps every normalized biallelic record passing QUAL, HQ depth
  (FORMAT/DP), and allele frequency (FORMAT/AD[0:1] / FORMAT/DP >= min_freq).
- Consensus VCF admits ONLY alleles appropriate for consensus assembly:
  * Haploid mode (--use-iupac false): VAF >= min_consensus_freq, GT=1.
  * IUPAC mode (--use-iupac true):
      - Indels/MNPs: VAF >= min_consensus_freq, GT=1/1 (never 0/1 IUPAC).
      - Biallelic SNPs: VAF >= max_iupac_freq -> GT=1/1 (fixed);
        min_iupac_freq <= VAF < max_iupac_freq -> GT=0/1 (IUPAC ambiguity).
  * Resolves same-position split-multiallelic collisions for bcftools consensus.
"""

import argparse
from collections import defaultdict
import sys


REQUIRED_INFO_HEADERS = [
    '##INFO=<ID=VAF,Number=1,Type=Float,Description="Per-allele variant allele frequency from FORMAT/AD and FORMAT/DP">',
    '##INFO=<ID=SUBCLONAL,Number=0,Type=Flag,Description="Variant below consensus allele frequency threshold">',
    '##INFO=<ID=CONSENSUS,Number=0,Type=Flag,Description="Variant included in consensus sequence">',
]


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare decoupled reported and consensus VCFs.")
    parser.add_argument("--input", required=True, help="Normalized input VCF path")
    parser.add_argument("--reported-out", required=True, help="Output path for reported VCF")
    parser.add_argument("--consensus-out", required=True, help="Output path for consensus VCF")
    parser.add_argument("--min-freq", type=float, required=True, help="Minimum VAF to report variant")
    parser.add_argument("--min-consensus-freq", type=float, default=0.50, help="Minimum VAF for haploid consensus / major indels")
    parser.add_argument("--min-coverage", type=int, required=True, help="Minimum HQ depth (FORMAT/DP)")
    parser.add_argument("--min-variant-qual", type=float, required=True, help="Minimum QUAL")
    parser.add_argument("--use-iupac", action="store_true", help="Enable IUPAC ambiguity genotypes for mixed biallelic SNPs")
    parser.add_argument("--min-iupac-freq", type=float, default=0.25, help="Lower VAF bound (inclusive) for IUPAC SNPs")
    parser.add_argument("--max-iupac-freq", type=float, default=0.75, help="Upper VAF bound (exclusive) for IUPAC SNPs")
    return parser.parse_args()


def extract_hq_depth_and_ad(parts):
    """Return (hq_dp, ref_ad, alt_ad) from FORMAT/DP+AD, falling back to INFO/DP4."""
    info_str = parts[7]
    info = {}
    for item in info_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            info[k] = v

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
                    ref_ad = ad_vals[0]
                    alt_ad = ad_vals[1]
                    if dp_str and dp_str != ".":
                        hq_dp = int(dp_str)
                    else:
                        hq_dp = sum(ad_vals)
                    return hq_dp, ref_ad, alt_ad
            except ValueError:
                pass

    if "DP4" in info:
        try:
            dp4 = [int(x) for x in info["DP4"].split(",")]
            if len(dp4) == 4:
                ref_ad = dp4[0] + dp4[1]
                alt_ad = dp4[2] + dp4[3]
                return ref_ad + alt_ad, ref_ad, alt_ad
        except ValueError:
            pass

    return 0, 0, 0


def set_sample_gt(parts, gt):
    """Set or prepend FORMAT/GT in VCF record parts."""
    if len(parts) < 10:
        return parts + ["GT", gt]
    fmt_keys = parts[8].split(":")
    smp_vals = parts[9].split(":")
    if "GT" in fmt_keys:
        idx = fmt_keys.index("GT")
        while len(smp_vals) <= idx:
            smp_vals.append(".")
        smp_vals[idx] = gt
    else:
        fmt_keys = ["GT"] + fmt_keys
        smp_vals = [gt] + smp_vals
    out = list(parts)
    out[8] = ":".join(fmt_keys)
    out[9] = ":".join(smp_vals)
    return out


def append_info_tags(info_str, vaf, is_consensus):
    """Append VAF and CONSENSUS/SUBCLONAL tags cleanly to INFO column."""
    items = [x for x in info_str.split(";") if x and x != "." and not x.startswith("VAF=") and x not in ("SUBCLONAL", "CONSENSUS")]
    items.append(f"VAF={vaf:.4f}")
    items.append("CONSENSUS" if is_consensus else "SUBCLONAL")
    return ";".join(items)


def prepare_vcfs(args):
    meta_lines = []
    chrom_header = None
    records = []

    with open(args.input, "r") as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith("##"):
                meta_lines.append(line)
            elif line.startswith("#CHROM"):
                chrom_header = line
            elif line.strip():
                records.append(line.split("\t"))

    if chrom_header is None:
        raise ValueError(f"Input VCF {args.input} is missing #CHROM header line.")

    existing_meta = "\n".join(meta_lines)
    for h_line in REQUIRED_INFO_HEADERS:
        tag_id = h_line.split("ID=")[1].split(",")[0]
        if f"##INFO=<ID={tag_id}," not in existing_meta:
            meta_lines.append(h_line)

    reported_lines = []
    consensus_candidates = defaultdict(list)
    pos_order = []

    for parts in records:
        if len(parts) < 8:
            continue
        chrom = parts[0]
        pos = int(parts[1])
        ref = parts[3]
        alt = parts[4]
        if alt == "." or "," in alt:
            # Skip non-variant or unexpected unsplit multiallelic records safely
            continue

        try:
            qual = float(parts[5]) if parts[5] != "." else 0.0
        except ValueError:
            qual = 0.0

        hq_dp, ref_ad, alt_ad = extract_hq_depth_and_ad(parts)
        if hq_dp <= 0:
            continue
        vaf = alt_ad / hq_dp

        if qual < args.min_variant_qual or hq_dp < args.min_coverage or vaf < args.min_freq:
            continue

        is_snp = (len(ref) == 1 and len(alt) == 1)

        if not args.use_iupac:
            eligible = vaf >= args.min_consensus_freq
            cons_gt = "1"
            rep_gt = "1" if eligible else "0/1"
        else:
            if not is_snp:
                eligible = vaf >= args.min_consensus_freq
                cons_gt = "1/1"
                rep_gt = "1/1" if eligible else "0/1"
            else:
                if vaf >= args.max_iupac_freq:
                    eligible = True
                    cons_gt = "1/1"
                elif vaf >= args.min_iupac_freq:
                    eligible = True
                    cons_gt = "0/1"
                else:
                    eligible = False
                    cons_gt = "0/1"
                rep_gt = cons_gt

        rep_parts = set_sample_gt(parts, rep_gt)
        rep_parts[7] = append_info_tags(rep_parts[7], vaf, eligible)
        reported_lines.append("\t".join(rep_parts))

        if eligible:
            key = (chrom, pos)
            if key not in consensus_candidates:
                pos_order.append(key)
            consensus_candidates[key].append({
                "parts": set_sample_gt(parts, cons_gt),
                "ref": ref,
                "alt": alt,
                "is_snp": is_snp,
                "vaf": vaf,
                "qual": qual,
                "hq_dp": hq_dp,
                "ref_ad": ref_ad,
                "alt_ad": alt_ad,
                "gt": cons_gt,
            })

    consensus_lines = []
    for key in pos_order:
        cands = consensus_candidates[key]
        if len(cands) == 1:
            c = cands[0]
            p = list(c["parts"])
            p[7] = append_info_tags(p[7], c["vaf"], True)
            consensus_lines.append("\t".join(p))
            continue

        # Multiple consensus-eligible alleles at the same (chrom, pos) after norm -m -any
        if (
            args.use_iupac
            and len(cands) == 2
            and all(c["is_snp"] and c["gt"] == "0/1" for c in cands)
            and (cands[0]["ref_ad"] / cands[0]["hq_dp"]) < args.min_iupac_freq
        ):
            # Both ALT SNPs are mixed while REF is minor -> emit biallelic-alt 1/2 for IUPAC code
            cands_sorted = sorted(cands, key=lambda x: (-x["vaf"], -x["qual"]))
            p = list(cands_sorted[0]["parts"])
            p[4] = f"{cands_sorted[0]['alt']},{cands_sorted[1]['alt']}"
            p = set_sample_gt(p, "1/2")
            # Adjust FORMAT/AD if present to ref,alt1,alt2
            if len(p) >= 10 and "AD" in p[8].split(":"):
                fmt_keys = p[8].split(":")
                smp_vals = p[9].split(":")
                ad_idx = fmt_keys.index("AD")
                smp_vals[ad_idx] = f"{cands_sorted[0]['ref_ad']},{cands_sorted[0]['alt_ad']},{cands_sorted[1]['alt_ad']}"
                p[9] = ":".join(smp_vals)
            p[7] = append_info_tags(p[7], cands_sorted[0]["vaf"] + cands_sorted[1]["vaf"], True)
            consensus_lines.append("\t".join(p))
        else:
            best = max(cands, key=lambda x: (x["vaf"], x["qual"]))
            p = list(best["parts"])
            p[7] = append_info_tags(p[7], best["vaf"], True)
            consensus_lines.append("\t".join(p))

    header_block = meta_lines + [chrom_header]
    with open(args.reported_out, "w") as f:
        f.write("\n".join(header_block + reported_lines) + "\n")
    with open(args.consensus_out, "w") as f:
        f.write("\n".join(header_block + consensus_lines) + "\n")


def main():
    args = parse_args()
    if not (0.0 <= args.min_freq <= 1.0):
        print("Error: --min-freq must be in [0, 1].", file=sys.stderr)
        sys.exit(1)
    if not (0.0 <= args.min_iupac_freq < args.max_iupac_freq <= 1.0):
        print("Error: require 0 <= --min-iupac-freq < --max-iupac-freq <= 1.", file=sys.stderr)
        sys.exit(1)
    prepare_vcfs(args)


if __name__ == "__main__":
    main()
