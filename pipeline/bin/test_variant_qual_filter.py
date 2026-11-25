#!/usr/bin/env python3
"""
Behavioral tests for variant QUAL/FORMAT-DP/VAF filtering, consensus VCF
decoupling, live bcftools filter execution, and nextflow_schema.json parity.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

BIN_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(BIN_DIR, "..", ".."))
sys.path.insert(0, BIN_DIR)

import prepare_consensus_vcf as pcv  # noqa: E402


EDGE_VCF = """##fileformat=VCFv4.2
##contig=<ID=SYN_VIRUS_1,length=2000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="HQ Depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSampleEdge
SYN_VIRUS_1\t10\t.\tA\tG\t19.9\tPASS\tDP=200\tGT:DP:AD\t1/1:200:10,190
SYN_VIRUS_1\t20\t.\tC\tT\t.\tPASS\tDP=200\tGT:DP:AD\t1/1:200:5,195
SYN_VIRUS_1\t30\t.\tG\tA\t20.0\tPASS\tDP=50\tGT:DP:AD\t1/1:50:1,49
SYN_VIRUS_1\t40\t.\tT\tC\t80.0\tPASS\tDP=200\tGT:DP:AD\t1/1:9:0,9
SYN_VIRUS_1\t50\t.\tA\tC\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:80,20
SYN_VIRUS_1\t60\t.\tC\tG\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:81,19
SYN_VIRUS_1\t70\t.\tG\tT\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:75,25
SYN_VIRUS_1\t80\t.\tT\tA\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:26,74
SYN_VIRUS_1\t90\t.\tA\tT\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:25,75
SYN_VIRUS_1\t100\t.\tG\tA\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:10,45
SYN_VIRUS_1\t100\t.\tG\tC\t80.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:10,45
SYN_VIRUS_1\t110\t.\tC\tA,G\t90.0\tPASS\tDP=100\tGT:DP:AD\t0/1:100:10,45,45
SYN_VIRUS_1\t120\t.\tA\t.\t90.0\tPASS\tDP=100\tGT:DP:AD\t0/0:100:100,0
"""


def parse_vcf_rows(path):
    rows = []
    with open(path, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                rows.append(line.rstrip("\n").split("\t"))
    return rows


class TestVariantQualAndSchema(unittest.TestCase):
    def test_qual_dp_af_boundaries_and_iupac_merge(self):
        with tempfile.TemporaryDirectory() as td:
            in_vcf = os.path.join(td, "in.vcf")
            rep_vcf = os.path.join(td, "rep.vcf")
            con_vcf = os.path.join(td, "con.vcf")
            with open(in_vcf, "w") as f:
                f.write(EDGE_VCF)

            # Haploid mode checks
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
            rep = {int(r[1]): r for r in parse_vcf_rows(rep_vcf) if int(r[1]) != 100}
            con = {int(r[1]): r for r in parse_vcf_rows(con_vcf)}

            # Rejected: QUAL 19.9 (pos 10), QUAL '.' (pos 20), FORMAT/DP=9 despite INFO/DP=200 (pos 40),
            # AF=0.19 < 0.20 (pos 60), unsplit multiallelic (pos 110), non-variant '.' (pos 120)
            for rejected in (10, 20, 40, 60, 110, 120):
                self.assertNotIn(rejected, rep)
                self.assertNotIn(rejected, con)

            # Accepted exact QUAL == 20.0 (pos 30) and exact min_freq == 0.20 (pos 50, subclonal)
            self.assertIn(30, rep)
            self.assertIn(30, con)
            self.assertIn(50, rep)
            self.assertIn("SUBCLONAL", rep[50][7].split(";"))
            self.assertNotIn(50, con)

            # IUPAC mode checks exact [0.25, 0.75) half-open bounds and 1/2 multi-alt merge at pos 100
            args_iupac = pcv.argparse.Namespace(**{**vars(args), "use_iupac": True})
            pcv.prepare_vcfs(args_iupac)
            con_iupac = {int(r[1]): r for r in parse_vcf_rows(con_vcf)}

            # pos 50 (0.20 < 0.25): excluded from consensus even in IUPAC mode
            self.assertNotIn(50, con_iupac)
            # pos 70 (0.25 inclusive): GT=0/1
            self.assertEqual(con_iupac[70][9].split(":")[0], "0/1")
            # pos 80 (0.74 < 0.75): GT=0/1
            self.assertEqual(con_iupac[80][9].split(":")[0], "0/1")
            # pos 90 (0.75 exclusive upper bound): GT=1/1 fixed!
            self.assertEqual(con_iupac[90][9].split(":")[0], "1/1")
            # pos 100 (REF 0.10 < 0.25, two 0.45 SNPs A & C): merged 1/2 A,C with 3-allele AD
            self.assertEqual(con_iupac[100][4], "A,C")
            self.assertEqual(con_iupac[100][9].split(":")[0], "1/2")
            self.assertEqual(con_iupac[100][9].split(":")[2], "10,45,45")

    def test_live_bcftools_filter_decoupled_from_base_qual(self):
        bcftools_bin = shutil.which("bcftools")
        if not bcftools_bin:
            self.skipTest("bcftools binary not available on PATH")

        with tempfile.TemporaryDirectory() as td:
            in_vcf = os.path.join(td, "raw.vcf")
            out_vcf = os.path.join(td, "filtered.vcf")
            with open(in_vcf, "w") as f:
                f.write(
                    "##fileformat=VCFv4.2\n"
                    "##contig=<ID=SYN_VIRUS_1,length=2000>\n"
                    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw DP\">\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"HQ DP\">\n"
                    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                    "SYN_VIRUS_1\t10\t.\tA\tG\t15.0\t.\tDP=100\tGT:DP:AD\t1/1:50:1,49\n"
                    "SYN_VIRUS_1\t20\t.\tC\tT\t25.0\t.\tDP=100\tGT:DP:AD\t1/1:8:0,8\n"
                    "SYN_VIRUS_1\t30\t.\tG\tA\t25.0\t.\tDP=100\tGT:DP:AD\t0/1:100:85,15\n"
                    "SYN_VIRUS_1\t40\t.\tT\tC\t25.0\t.\tDP=100\tGT:DP:AD\t0/1:100:70,30\n"
                )

            expr = "QUAL>=20 && FORMAT/DP>=10 && FORMAT/DP>0 && (FORMAT/AD[0:1]/FORMAT/DP)>=0.2"
            subprocess.run(
                [bcftools_bin, "filter", "-i", expr, "-Ov", "-o", out_vcf, in_vcf],
                check=True,
                capture_output=True,
                text=True,
            )
            kept = [int(r[1]) for r in parse_vcf_rows(out_vcf)]
            self.assertEqual(kept, [40], f"Expected only pos 40 to pass live bcftools filter, got {kept}")

    def test_nextflow_schema_and_config_parity(self):
        schema_path = os.path.join(REPO_ROOT, "nextflow_schema.json")
        config_path = os.path.join(REPO_ROOT, "pipeline", "nextflow.config")
        self.assertTrue(os.path.exists(schema_path), f"Missing {schema_path}")

        with open(schema_path, "r") as f:
            schema = json.load(f)

        schema_params = set()
        for sec in schema.get("$defs", {}).values():
            schema_params.update(sec.get("properties", {}).keys())

        with open(config_path, "r") as f:
            cfg_text = f.read()

        # Extract params { ... } block keys
        in_params = False
        brace_depth = 0
        config_params = set()
        for line in cfg_text.splitlines():
            stripped = line.strip()
            if not in_params:
                if re.match(r"^params\s*\{", stripped):
                    in_params = True
                    brace_depth = 1
                continue
            brace_depth += stripped.count("{") - stripped.count("}")
            if brace_depth <= 0:
                break
            m = re.match(r"^([a-zA-Z_][a-zA-Z0-9_]*)\s*=", stripped)
            if m:
                config_params.add(m.group(1))

        self.assertGreater(len(config_params), 20)
        self.assertEqual(
            config_params,
            schema_params,
            f"Config vs nextflow_schema mismatch: only_in_config={config_params - schema_params}, "
            f"only_in_schema={schema_params - config_params}",
        )


if __name__ == "__main__":
    unittest.main()
