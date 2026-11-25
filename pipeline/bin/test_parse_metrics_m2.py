#!/usr/bin/env python3
"""
Milestone 2 unit test suite covering:
- parse_fastp (Illumina fastp + ONT fastplong, null duplication_rate, fail-loudly)
- parse_vcf (FORMAT/AD+DP priority, FILTER enforcement, missing QUAL/DP/AF non-fabrication)
- Two-pass primer pairing (ARTIC _LEFT_1/-alt1, pre-direction RSV_S1_01_LEFT+RSV_S1_RIGHT, Pan-CoV-F/R)
- Callable-target QC tiers & and-gated amplicon pass rule
- Exact contig lookup in scripts/retrieve_hits.py
- JSON Schema validation against schema/viromics_run_report.schema.json
"""

import json
import os
import sys
import tempfile
import unittest

BIN_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(BIN_DIR, "..", ".."))
sys.path.insert(0, BIN_DIR)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))

import parse_metrics  # noqa: E402
import retrieve_hits  # noqa: E402


class TestMilestone2Engine(unittest.TestCase):
    def test_parse_fastp_illumina_and_fastplong_ont(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            illumina_json = os.path.join(tmpdir, "s1_fastp.json")
            with open(illumina_json, "w") as f:
                json.dump({
                    "summary": {
                        "before_filtering": {"total_reads": 1000, "total_bases": 150000, "q30_rate": 0.91, "read1_mean_length": 150},
                        "after_filtering": {"total_reads": 950, "total_bases": 141000, "q30_rate": 0.94, "read1_mean_length": 148},
                    },
                    "duplication": {"rate": 0.045},
                }, f)

            qc_ill = parse_metrics.parse_fastp(illumina_json)
            self.assertEqual(qc_ill["raw_reads"], 1000)
            self.assertEqual(qc_ill["duplication_rate"], 0.045)
            self.assertEqual(qc_ill["q30_rate_trimmed"], 0.94)
            self.assertEqual(qc_ill["mean_read_length_raw"], 150.0)

            ont_json = os.path.join(tmpdir, "b1_fastp.json")
            with open(ont_json, "w") as f:
                json.dump({
                    "summary": {
                        "before_filtering": {"total_reads": 5000, "total_bases": 2100000, "q30_rate": 0.55, "read_mean_length": 420},
                        "after_filtering": {"total_reads": 4500, "total_bases": 1900000, "q30_rate": 0.61, "read_mean_length": 422},
                    }
                }, f)

            qc_ont = parse_metrics.parse_fastp(ont_json)
            self.assertIsNone(qc_ont["duplication_rate"])
            self.assertEqual(qc_ont["q30_rate_raw"], 0.55)
            self.assertEqual(qc_ont["mean_read_length_trimmed"], 422.0)

            with self.assertRaises(FileNotFoundError):
                parse_metrics.parse_fastp(os.path.join(tmpdir, "missing.json"))

    def test_parse_vcf_format_ad_dp_filter_and_no_fabricated_defaults(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf_path = os.path.join(tmpdir, "test.vcf")
            with open(vcf_path, "w") as f:
                f.writelines([
                    "##fileformat=VCFv4.2\n",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
                    # Soft-filtered record must be excluded
                    "ref1\t10\t.\tA\tG\t50\tLowQual\tDP=200\tGT:DP:AD\t0/1:200:50,150\n",
                    # Missing QUAL ('.') stays None; FORMAT/AD=80,20 DP=100 -> AF=0.20 (even if INFO/DP4 differed)
                    "ref1\t20\t.\tC\tT\t.\tPASS\tDP=999;DP4=0,0,0,0;SUBCLONAL\tGT:DP:AD\t0/1:100:80,20\n",
                    # Missing DP/AF in INFO and no FORMAT -> dp=0, filtered out instead of fabricating dp=100, af=1.0
                    "ref1\t30\t.\tG\tA\t60\tPASS\t.\n",
                ])

            variants = parse_metrics.parse_vcf(vcf_path, {}, ref_seqs={"ref1": "A" * 100}, min_freq=0.2)
            self.assertEqual(len(variants), 1)
            v = variants[0]
            self.assertEqual(v["position"], 20)
            self.assertIsNone(v["quality"])
            self.assertEqual(v["depth"], 100)
            self.assertEqual(v["frequency"], 0.2)
            self.assertFalse(v["is_consensus"])

    def test_primer_pairing_variants_and_and_gated_amplicons(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = os.path.join(tmpdir, "primers.bed")
            with open(bed_path, "w") as f:
                f.writelines([
                    # ARTIC alt index after direction
                    "SARS2\t30\t54\tnCoV-2019_1_LEFT_1\t1\t+\n",
                    "SARS2\t35\t58\tnCoV-2019_1_LEFT_alt1\t1\t+\n",
                    "SARS2\t400\t424\tnCoV-2019_1_RIGHT\t1\t-\n",
                    # Pre-direction alt index (real RSV-Eden style)
                    "RSV\t10\t32\tRSV_S1_01_LEFT\t1\t+\n",
                    "RSV\t500\t522\tRSV_S1_RIGHT\t1\t-\n",
                    # Pan-CoV F/R style
                    "FIPV\t0\t21\tPan-CoV-F\t1\t+\n",
                    "FIPV\t414\t434\tPan-CoV-R\t1\t-\n",
                ])

            amps = {a["chrom"]: a for a in parse_metrics.parse_primer_amplicons(bed_path)}
            # SARS2 insert uses max left end (58) to min right start (400)
            self.assertEqual((amps["SARS2"]["start"], amps["SARS2"]["end"]), (58, 400))
            self.assertEqual((amps["RSV"]["start"], amps["RSV"]["end"]), (32, 500))
            self.assertEqual((amps["FIPV"]["start"], amps["FIPV"]["end"]), (21, 414))

            # Amplicon with 50% complete dropout (0x) and 50% at 1000x has mean_depth=500 >= 50,
            # pct_above=50.0 < 90.0 -> MUST FAIL under and-gated rule!
            depths = {"FIPV": [0] * 21 + [0] * 197 + [1000] * 196 + [0] * 20}
            amp_stats = parse_metrics.calculate_amplicon_coverage([amps["FIPV"]], depths, min_cov=50)
            self.assertEqual(len(amp_stats), 1)
            self.assertEqual(amp_stats[0]["status"], "fail")

    def test_callable_target_qc_tiers_pancov(self):
        # 434-bp FIPV amplicon: terminal 21 + 20 = 41 bp N-masked after primer clip, 393 bp insert at 200x -> PASS
        ref_info = [{"name": "FIPV_NC_002306", "length": 434, "sequence": "A" * 434}]
        primer_amps = [{"chrom": "FIPV_NC_002306", "start": 21, "end": 414, "name": "Pan-CoV"}]
        qc_thresholds = {
            "pass_depth": 100, "pass_breadth": 95.0, "pass_max_n": 5.0,
            "warn_depth": 25, "warn_breadth": 80.0, "warn_max_n": 10.0,
        }
        depths_pass = [0] * 21 + [200] * 393 + [0] * 20
        s_pass = {
            "sample_id": "barcode01",
            "full_depths": {"FIPV_NC_002306": depths_pass},
            "coverage": {
                "FIPV_NC_002306": {
                    "mean_depth": round(sum(depths_pass) / 434, 2),
                    "coverage_breadth_10x": 90.55,
                    "coverage_breadth_30x": 90.55,
                    "coverage_breadth_min_cov": 90.55,
                    "mapped_reads": 500,
                    "coverage_profile": [],
                    "low_coverage_intervals": [[1, 21], [415, 434]],
                }
            },
            "consensus": {"FIPV_NC_002306": {"length": 434, "gc_content": 45.0, "n_content": 9.45, "sequence": "N" * 21 + "A" * 393 + "N" * 20}},
        }
        parse_metrics.annotate_callable_targets_and_qc(
            s_pass, ref_info, primer_amps, [], min_cov=50, mask_depth=50, qc_thresholds=qc_thresholds
        )
        self.assertEqual(s_pass["qc_status"], "PASS")
        self.assertEqual(s_pass["coverage"]["FIPV_NC_002306"]["callable_target_length"], 393)
        self.assertEqual(s_pass["coverage"]["FIPV_NC_002306"]["callable_breadth_min_cov"], 100.0)
        self.assertEqual(s_pass["consensus"]["FIPV_NC_002306"]["internal_n_content"], 0.0)
        self.assertEqual(s_pass["consensus"]["FIPV_NC_002306"]["terminal_n_count"], 41)

    def test_retrieve_hits_exact_contig_matching(self):
        records = {"S1_CoV_10": "GGGG", "S1_CoV_1": "AAAA"}
        self.assertEqual(retrieve_hits.find_contig_sequence(records, "S1", "CoV_1"), "AAAA")
        # If CoV_1 is absent, must return None rather than matching CoV_10 by substring!
        self.assertIsNone(retrieve_hits.find_contig_sequence({"S1_CoV_10": "GGGG"}, "S1", "CoV_1"))

    def test_mock_report_matches_json_schema_file(self):
        mock_path = os.path.join(REPO_ROOT, "dashboard", "public", "mock_run_report.json")
        schema_path = os.path.join(REPO_ROOT, "schema", "viromics_run_report.schema.json")
        with open(mock_path, "r") as f:
            report = json.load(f)
        with open(schema_path, "r") as f:
            schema = json.load(f)

        parse_metrics.validate_report_schema(report)
        try:
            import jsonschema
            jsonschema.validate(instance=report, schema=schema)
        except ImportError:
            self.skipTest("jsonschema package not installed in this Python interpreter")


if __name__ == "__main__":
    unittest.main()
