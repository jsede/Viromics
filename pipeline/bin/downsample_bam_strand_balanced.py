#!/usr/bin/env python3
"""
Deterministic, strand-balanced BAM downsampler for ultra-deep viral amplicons.

Prevents coordinate-sorted pileup depth caps (bcftools mpileup -d) from
sampling 100% forward reads and 0 reverse reads when tens of thousands of
amplicon reads share identical clipped start positions.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile


def parse_args():
    parser = argparse.ArgumentParser(description="Deterministic strand-balanced BAM downsampler.")
    parser.add_argument("--bam", required=True, help="Input coordinate-sorted, indexed BAM")
    parser.add_argument("--out", required=True, help="Output BAM path")
    parser.add_argument("--max-depth", type=int, default=100000, help="Max total depth cap per contig (<=0 disables)")
    parser.add_argument("--seed", type=int, default=42, help="Deterministic RNG seed for samtools view --subsample-seed")
    parser.add_argument("--map-qual", type=int, default=20, help="MAPQ threshold matching mpileup -q")
    parser.add_argument("--base-qual", type=int, default=15, help="BaseQ threshold matching mpileup -Q")
    parser.add_argument("--threads", type=int, default=1, help="Samtools threads")
    return parser.parse_args()


def max_depth_per_contig(bam_path, map_qual, base_qual, contig=None):
    cmd = [
        "samtools", "depth", "-J",
        "-Q", str(map_qual),
        "-q", str(base_qual),
        "-G", "0xF04",
    ]
    if contig:
        cmd.extend(["-r", contig])
    cmd.append(bam_path)
    res = subprocess.run(cmd, capture_output=True, text=True, check=True)
    max_by_contig = {}
    for line in res.stdout.splitlines():
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        c = parts[0]
        d = int(parts[2])
        if d > max_by_contig.get(c, 0):
            max_by_contig[c] = d
    return max_by_contig


def count_bam_records(bam_path):
    res = subprocess.run(["samtools", "view", "-c", bam_path], capture_output=True, text=True, check=True)
    return int(res.stdout.strip())


def copy_and_index(src_bam, dst_bam, threads):
    shutil.copy2(src_bam, dst_bam)
    subprocess.run(["samtools", "index", "-@", str(threads), dst_bam], check=True)


def main():
    args = parse_args()
    if args.max_depth <= 0:
        copy_and_index(args.bam, args.out, args.threads)
        return

    max_by_contig = max_depth_per_contig(args.bam, args.map_qual, args.base_qual)
    if not max_by_contig or max(max_by_contig.values()) <= args.max_depth:
        copy_and_index(args.bam, args.out, args.threads)
        return

    idx_res = subprocess.run(["samtools", "idxstats", args.bam], capture_output=True, text=True, check=True)
    contigs = []
    for line in idx_res.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 3 and parts[0] != "*" and int(parts[2]) > 0:
            contigs.append(parts[0])

    if not contigs:
        copy_and_index(args.bam, args.out, args.threads)
        return

    target_per_strand = args.max_depth / 2.0

    with tempfile.TemporaryDirectory() as tmpdir:
        pieces = []
        for idx, contig in enumerate(contigs):
            obs_max = max_by_contig.get(contig, 0)
            if obs_max <= args.max_depth:
                piece = os.path.join(tmpdir, f"c{idx}_full.bam")
                subprocess.run(
                    ["samtools", "view", "-@", str(args.threads), "-b", "-F", "0x904", args.bam, contig, "-o", piece],
                    check=True,
                )
                if count_bam_records(piece) > 0:
                    pieces.append(piece)
                continue

            for strand_label, flag_args in [("fwd", ["-F", "0x914"]), ("rev", ["-f", "0x10", "-F", "0x904"])]:
                raw_strand = os.path.join(tmpdir, f"c{idx}_{strand_label}.bam")
                subprocess.run(
                    ["samtools", "view", "-@", str(args.threads), "-b"] + flag_args + [args.bam, contig, "-o", raw_strand],
                    check=True,
                )
                if count_bam_records(raw_strand) == 0:
                    continue
                subprocess.run(["samtools", "index", raw_strand], check=True)
                s_depths = max_depth_per_contig(raw_strand, args.map_qual, args.base_qual, contig=contig)
                strand_max = s_depths.get(contig, 0)
                if strand_max > target_per_strand:
                    frac = max(0.0001, min(0.999999, target_per_strand / strand_max))
                    sub_strand = os.path.join(tmpdir, f"c{idx}_{strand_label}_sub.bam")
                    subprocess.run(
                        [
                            "samtools", "view", "-@", str(args.threads), "-b",
                            "--subsample-seed", str(args.seed),
                            "--subsample", f"{frac:.6f}",
                            raw_strand, "-o", sub_strand,
                        ],
                        check=True,
                    )
                    print(
                        f"[downsample] {contig} ({strand_label}): max_depth={strand_max} > cap={target_per_strand:.0f}; "
                        f"subsampled frac={frac:.5f} (seed={args.seed})",
                        file=sys.stderr,
                    )
                    if count_bam_records(sub_strand) > 0:
                        pieces.append(sub_strand)
                else:
                    pieces.append(raw_strand)

        if not pieces:
            copy_and_index(args.bam, args.out, args.threads)
            return

        if len(pieces) == 1:
            subprocess.run(["samtools", "sort", "-@", str(args.threads), "-o", args.out, pieces[0]], check=True)
        else:
            merged = os.path.join(tmpdir, "merged.bam")
            subprocess.run(["samtools", "merge", "-f", "-@", str(args.threads), merged] + pieces, check=True)
            subprocess.run(["samtools", "sort", "-@", str(args.threads), "-o", args.out, merged], check=True)
        subprocess.run(["samtools", "index", "-@", str(args.threads), args.out], check=True)


if __name__ == "__main__":
    main()
