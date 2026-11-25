#!/usr/bin/env python3
"""
Deterministic synthetic viral amplicon WGS fixture generator for CI & local E2E testing.

Generates a compact (<100 KB total) fixture under pipeline/test_data/ with planted ground truth:
  - SYN_VIRUS_1 (1,430 bp) + ORF1 (101..1390, + strand, phase=0)
  - 6 tiled 25-bp primer pairs (Synth_1..Synth_6, 280-bp amplicons, contiguous inserts)
  - Samples (Illumina PE 146-150 bp & ONT 280 bp):
      * SynthPos:     105x+ coverage (PASS), clonal SNP 200 G>A (ORF1: A34T),
                      ~31% subclonal 1-bp insertion 370 T>TA (reported in Illumina VCF,
                      strictly excluded from consensus in BOTH Illumina & ONT),
                      clonal 6-bp in-frame deletion 572..577 (pos 571 CATGCTG>C, consensus 1424 bp),
                      and decoy mismatch at pos 10 inside Synth_1_LEFT (must be clipped cleanly).
      * SynthDropout: Drops Amp 4 (230-bp dropout gap 716..945 -> N-masked & low_coverage_intervals),
                      WT at pos 200 (1 SNP distance from SynthPos), shares clonal 6-bp deletion.
      * SynthNTC:     Empty 0-read FASTQ control (completes without crashing, classifies FAIL).
"""

import gzip
import os
import random

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
SEED = 2026


def comp(b):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(b, "N")


def revcomp(seq):
    return "".join(comp(b) for b in reversed(seq))


def build_reference():
    rng = random.Random(SEED)
    bases = ["A", "C", "G", "T"]
    seq = []
    for _ in range(1430):
        choices = [b for b in bases if not (len(seq) >= 2 and seq[-1] == b and seq[-2] == b)]
        seq.append(rng.choice(choices))

    # ORF1 runs 1-based 101..1390 (0-based 100..1390, 430 codons)
    seq[100:103] = list("ATG")
    seq[1387:1390] = list("TAA")

    # Eliminate internal stop codons in ORF1 codons 2..429
    for idx0 in range(103, 1387, 3):
        codon = "".join(seq[idx0 : idx0 + 3])
        if codon in ("TAA", "TAG", "TGA"):
            seq[idx0 + 2] = "C"

    # Planted clonal SNP locus: codon 34 at 1-based 200..202 (0-based 199..202) = GCT (Ala)
    # 200 G->A yields ACT (Thr) -> A34T
    seq[196:205] = list("CGAGCTGAT")

    # Planted ~31% subclonal insertion locus: 1-based pos 370 (0-based 369) = 'T', flanked by C/G
    # Inside Amp 2 R1+R2 overlap [360..380]; inserting 'A' after 370 gives CGTACGA -> 370 T>TA
    seq[366:374] = list("AGCTCGAT")

    # Planted clonal 6-bp in-frame deletion locus: 1-based 572..577 (0-based 571..577)
    # Codon 157 (569..571, 0-based 568..571) = GAC; Codons 158-159 (572..577, 0-based 571..577) = ATGCTG
    # Deleting 572..577 leaves anchor 571 C -> pos 571 REF=CATGCTG ALT=C
    seq[568:580] = list("GACATGCTGATC")

    return "".join(seq)


AMPLICONS = [
    ("Synth_1", 0, 25, 255, 280, 1),
    ("Synth_2", 230, 255, 485, 510, 2),
    ("Synth_3", 460, 485, 715, 740, 1),
    ("Synth_4", 690, 715, 945, 970, 2),
    ("Synth_5", 920, 945, 1175, 1200, 1),
    ("Synth_6", 1150, 1175, 1405, 1430, 2),
]


def mutate_amplicon(ref_seq, start, end, clonal_snp=False, subclonal_ins=False, clonal_del=False, primer_decoy=False):
    """Build mutated slice matching [start, end) in reference coordinates."""
    out = []
    for idx0 in range(start, end):
        pos1 = idx0 + 1
        base = ref_seq[idx0]
        if primer_decoy and pos1 == 10:
            # Decoy mismatch inside Synth_1_LEFT (1-based 1..25)
            out.append("T" if base != "T" else "G")
        elif clonal_snp and pos1 == 200:
            out.append("A")
        elif clonal_del and 572 <= pos1 <= 577:
            continue
        elif subclonal_ins and pos1 == 370:
            out.append(base)
            out.append("A")
        else:
            out.append(base)
    return "".join(out)


def write_fastq_gz(path, records):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt", compresslevel=9) as gz:
        for hdr, seq, qual in records:
            gz.write(f"@{hdr}\n{seq}\n+\n{qual}\n")


def generate_illumina_sample(ref_seq, amp_indices, clonal_snp=False, subclonal_ins=False, clonal_del=False):
    r1_recs = []
    r2_recs = []
    total_reads = 105
    sub_count = 33  # 33/105 = 31.43% subclonal fraction (~31%)

    counter = 0
    for amp_idx in amp_indices:
        amp_name, l_start, _l_end, _r_start, r_end, _pool = AMPLICONS[amp_idx]
        for i in range(total_reads):
            use_ins = subclonal_ins and (amp_name == "Synth_2") and (i < sub_count)
            amp_seq = mutate_amplicon(
                ref_seq,
                l_start,
                r_end,
                clonal_snp=clonal_snp,
                subclonal_ins=use_ins,
                clonal_del=clonal_del,
                primer_decoy=(amp_name == "Synth_1"),
            )
            # Mild 3'-end length jitter (147..150 bp) creates realistic VDB variance around overlap midpoint
            r1_len = 150 - (i % 4)
            r2_len = 150 - ((i + 2) % 4)
            counter += 1
            r1 = amp_seq[:r1_len]
            r2 = revcomp(amp_seq[-r2_len:])
            r1_recs.append((f"r{counter}/1", r1, "F" * len(r1)))
            r2_recs.append((f"r{counter}/2", r2, "F" * len(r2)))

    return r1_recs, r2_recs


def generate_ont_sample(ref_seq, amp_indices, clonal_snp=False, subclonal_ins=False, clonal_del=False):
    recs = []
    total_reads = 105
    sub_count = 33  # 33/105 = 31.43% subclonal fraction (~31%)

    counter = 0
    for amp_idx in amp_indices:
        amp_name, l_start, _l_end, _r_start, r_end, _pool = AMPLICONS[amp_idx]
        for i in range(total_reads):
            use_ins = subclonal_ins and (amp_name == "Synth_2") and (i < sub_count)
            amp_seq = mutate_amplicon(
                ref_seq,
                l_start,
                r_end,
                clonal_snp=clonal_snp,
                subclonal_ins=use_ins,
                clonal_del=clonal_del,
                primer_decoy=(amp_name == "Synth_1"),
            )
            counter += 1
            seq_out = amp_seq if (i % 2 == 0) else revcomp(amp_seq)
            qual_out = "?" * len(seq_out)
            recs.append((f"r{counter}", seq_out, qual_out))

    return recs


def main():
    ref_seq = build_reference()

    fasta_path = os.path.join(OUT_DIR, "ref_synth.fasta")
    with open(fasta_path, "w") as f:
        f.write(">SYN_VIRUS_1 Synthetic viral genome for CI E2E testing\n")
        for i in range(0, len(ref_seq), 80):
            f.write(ref_seq[i : i + 80] + "\n")

    gff_path = os.path.join(OUT_DIR, "ref_synth.gff3")
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n")
        f.write("##sequence-region SYN_VIRUS_1 1 1430\n")
        f.write("SYN_VIRUS_1\tsynthetic\tgene\t101\t1390\t.\t+\t.\tID=gene-ORF1;Name=ORF1;gene=ORF1\n")
        f.write("SYN_VIRUS_1\tsynthetic\tCDS\t101\t1390\t.\t+\t0\tID=cds-ORF1;Parent=gene-ORF1;Name=ORF1;gene=ORF1;protein_id=YP_SYN001.1;product=synthetic polyprotein\n")

    bed_path = os.path.join(OUT_DIR, "primers_synth.bed")
    with open(bed_path, "w") as f:
        for amp_name, l_s, l_e, r_s, r_e, pool in AMPLICONS:
            f.write(f"SYN_VIRUS_1\t{l_s}\t{l_e}\t{amp_name}_LEFT\t{pool}\t+\n")
            f.write(f"SYN_VIRUS_1\t{r_s}\t{r_e}\t{amp_name}_RIGHT\t{pool}\t-\n")

    meta_path = os.path.join(OUT_DIR, "metadata_synth.csv")
    with open(meta_path, "w") as f:
        f.write("sample_id,collection_date,location,host\n")
        f.write("SynthPos,2026-05-01,Melbourne,Human\n")
        f.write("SynthDropout,2026-05-02,Sydney,Human\n")
        f.write("SynthNTC,2026-05-01,Melbourne,Control\n")

    all_amps = [0, 1, 2, 3, 4, 5]
    dropout_amps = [0, 1, 2, 4, 5]  # Omit Amp 4 -> dropout 716..945

    samples_spec = [
        ("SynthPos", all_amps, True, True, True),
        ("SynthDropout", dropout_amps, False, False, True),
        ("SynthNTC", [], False, False, False),
    ]

    ill_sheet_rows = ["sample_id,platform,fastq_1,fastq_2"]
    ont_sheet_rows = ["sample_id,platform,fastq_1"]

    for sid, amps, c_snp, s_ins, c_del in samples_spec:
        if amps:
            r1, r2 = generate_illumina_sample(ref_seq, amps, c_snp, s_ins, c_del)
            ont_r = generate_ont_sample(ref_seq, amps, c_snp, s_ins, c_del)
        else:
            r1, r2, ont_r = [], [], []

        write_fastq_gz(os.path.join(OUT_DIR, "illumina", f"{sid}_R1.fastq.gz"), r1)
        write_fastq_gz(os.path.join(OUT_DIR, "illumina", f"{sid}_R2.fastq.gz"), r2)
        write_fastq_gz(os.path.join(OUT_DIR, "ont", f"{sid}.fastq.gz"), ont_r)

        ill_sheet_rows.append(
            f"{sid},illumina,pipeline/test_data/illumina/{sid}_R1.fastq.gz,pipeline/test_data/illumina/{sid}_R2.fastq.gz"
        )
        ont_sheet_rows.append(f"{sid},ont,pipeline/test_data/ont/{sid}.fastq.gz")

    with open(os.path.join(OUT_DIR, "samplesheet_illumina.csv"), "w") as f:
        f.write("\n".join(ill_sheet_rows) + "\n")

    with open(os.path.join(OUT_DIR, "samplesheet_ont.csv"), "w") as f:
        f.write("\n".join(ont_sheet_rows) + "\n")

    print("Generated synthetic fixture under pipeline/test_data/.")


if __name__ == "__main__":
    main()
