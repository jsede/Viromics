#!/usr/bin/env python3
import os
import sys
import tempfile
import unittest

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import prepare_consensus_vcf as pcv


SAMPLE_VCF = """##fileformat=VCFv4.2
##contig=<ID=NC_045512.2,length=29870>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="HQ Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample03
NC_045512.2\t5673\t.\tC\tCT\t205.3\tPASS\tINDEL;DP=380\tGT:DP:AD\t0/1:371:254,117
NC_045512.2\t11287\t.\tGTCTGGTTTT\tG\t228.4\tPASS\tINDEL;DP=1260\tGT:DP:AD\t1/1:1250:52,1198
NC_045512.2\t11750\t.\tC\tT\t221.4\tPASS\tDP=1050\tGT:DP:AD\t0/1:1035:746,289
NC_045512.2\t14408\t.\tC\tT\t228.4\tPASS\tDP=900\tGT:DP:AD\t1/1:890:10,880
NC_045512.2\t20000\t.\tA\tG\t150.0\tPASS\tDP=200\tGT:DP:AD\t0/1:180:162,18
NC_045512.2\t25000\t.\tG\tA\t90.0\tPASS\tDP=80\tGT:DP:AD\t0/1:80:48,32
NC_045512.2\t25000\t.\tG\tC\t95.0\tPASS\tDP=80\tGT:DP:AD\t0/1:80:40,40
"""


def read_records(vcf_path):
    rows = []
    with open(vcf_path) as f:
        for l in f:
            if not l.startswith("#") and l.strip():
                rows.append(l.strip().split("\t"))
    return rows


class TestPrepareConsensusVCF(unittest.TestCase):
    def test_haploid_excludes_subclonal_indel_and_snp(self):
        with tempfile.TemporaryDirectory() as td:
            in_vcf = os.path.join(td, "in.vcf")
            rep_vcf = os.path.join(td, "rep.vcf")
            con_vcf = os.path.join(td, "con.vcf")
            with open(in_vcf, "w") as f:
                f.write(SAMPLE_VCF)

            args = pcv.argparse.Namespace(
                input=in_vcf,
                reported_out=rep_vcf,
                consensus_out=con_vcf,
                min_freq=0.20,
                min_consensus_freq=0.50,
                min_coverage=10,
                min_variant_qual=20.0,
                use_iupac=False,
                min_iupac_freq=0.25,
                max_iupac_freq=0.75,
            )
            pcv.prepare_vcfs(args)

            rep = read_records(rep_vcf)
            con = read_records(con_vcf)

            # pos 20000 has AF=18/180=0.10 < 0.20 -> excluded from both
            rep_positions = [(int(r[1]), r[4]) for r in rep]
            self.assertIn((5673, "CT"), rep_positions)
            self.assertIn((11750, "T"), rep_positions)
            self.assertNotIn((20000, "G"), rep_positions)

            # Consensus must NOT contain pos 5673 (AF=0.315) or pos 11750 (AF=0.279)
            con_positions = [(int(r[1]), r[4], r[9].split(":")[0]) for r in con]
            self.assertEqual(
                con_positions,
                [
                    (11287, "G", "1"),
                    (14408, "T", "1"),
                    (25000, "C", "1"),  # conflict resolved to highest VAF (0.50 vs 0.40)
                ],
            )

    def test_iupac_mode_gates_snps_vs_indels(self):
        with tempfile.TemporaryDirectory() as td:
            in_vcf = os.path.join(td, "in.vcf")
            rep_vcf = os.path.join(td, "rep.vcf")
            con_vcf = os.path.join(td, "con.vcf")
            with open(in_vcf, "w") as f:
                f.write(SAMPLE_VCF)

            args = pcv.argparse.Namespace(
                input=in_vcf,
                reported_out=rep_vcf,
                consensus_out=con_vcf,
                min_freq=0.20,
                min_consensus_freq=0.50,
                min_coverage=10,
                min_variant_qual=20.0,
                use_iupac=True,
                min_iupac_freq=0.25,
                max_iupac_freq=0.75,
            )
            pcv.prepare_vcfs(args)

            con = read_records(con_vcf)
            con_map = {int(r[1]): (r[4], r[9].split(":")[0]) for r in con}

            # Subclonal indel 5673 (AF=0.315) STILL excluded in IUPAC mode!
            self.assertNotIn(5673, con_map)
            # Major deletion 11287 (AF=0.958) -> GT=1/1
            self.assertEqual(con_map[11287], ("G", "1/1"))
            # Mixed SNP 11750 (AF=0.279 in [0.25, 0.75)) -> GT=0/1 (IUPAC Y)
            self.assertEqual(con_map[11750], ("T", "0/1"))
            # Fixed SNP 14408 (AF=0.989 >= 0.75) -> GT=1/1 (NOT 0/1)
            self.assertEqual(con_map[14408], ("T", "1/1"))
            # Conflict at 25000 with two 0/1 SNPs (0.50 + 0.40, ref_ad=40/80=0.50 >= 0.25) -> best single 0/1
            self.assertEqual(con_map[25000], ("C", "0/1"))


if __name__ == "__main__":
    unittest.main()
