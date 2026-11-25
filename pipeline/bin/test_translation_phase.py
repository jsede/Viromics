#!/usr/bin/env python3
"""
Biological ground-truth tests for GFF3 phase offsets, negative-strand codons,
partial-codon bounds, real SARS-CoV-2 ORF1a/ORF1ab ribosomal slippage,
-1/+1 curated frameshift logic, coding indel classification, and overlapping genes.
"""

import os
import sys
import tempfile
import unittest

BIN_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(BIN_DIR, "..", ".."))
sys.path.insert(0, BIN_DIR)
import parse_metrics  # noqa: E402


class TestTranslationPhase(unittest.TestCase):
    def test_positive_strand_phase_1_biological_ground_truth(self):
        # GFF3 spec: CDS 100..150 on + strand with phase=1 skips base 100 (1-based),
        # so Codon 1 = 101..103 (0-based 100..103), Codon 2 = 104..106 (0-based 103..106).
        ref_seq_list = ["N"] * 200
        ref_seq_list[99] = "G"   # pos 100: skipped partial codon base
        ref_seq_list[100] = "A"  # pos 101: Codon 1 base 1
        ref_seq_list[101] = "T"  # pos 102: Codon 1 base 2
        ref_seq_list[102] = "G"  # pos 103: Codon 1 base 3 (ATG = M)
        ref_seq_list[103] = "A"  # pos 104: Codon 2 base 1
        ref_seq_list[104] = "C"  # pos 105: Codon 2 base 2 (ACT = T)
        ref_seq_list[105] = "T"  # pos 106: Codon 2 base 3
        ref_seq = "".join(ref_seq_list)

        gene = {
            "name": "GeneP1",
            "type": "CDS",
            "start": 100,
            "end": 150,
            "strand": "+",
            "phase": 1,
            "multi_part": False,
        }

        # pos 100 is inside the skipped leading partial base -> fail closed
        partial_res = parse_metrics.get_aa_change(100, "G", "A", gene, ref_seq)
        self.assertEqual(partial_res["aa_change"], "AA annotation unavailable")
        self.assertEqual(partial_res["status"], "incomplete_partial_codon")

        # pos 103 (G->A) mutates ATG (M1) -> ATA (I1): M1I
        res_c1 = parse_metrics.get_aa_change(103, "G", "A", gene, ref_seq)
        self.assertEqual(res_c1["aa_change"], "M1I")
        self.assertEqual(res_c1["effect"], "missense")
        self.assertIsNone(res_c1["status"])

        # pos 105 (C->G) mutates ACT (T2) -> AGT (S2): T2S
        res_c2 = parse_metrics.get_aa_change(105, "C", "G", gene, ref_seq)
        self.assertEqual(res_c2["aa_change"], "T2S")
        self.assertEqual(res_c2["effect"], "missense")
        self.assertIsNone(res_c2["status"])

    def test_positive_strand_phase_2_and_trailing_partial_codon(self):
        # CDS 10..17 (len 8), phase=2 skips 10,11; Codon 1 = 12..14; Codon 2 = 15..17 (exact!)
        # If CDS is 10..16 (len 7), phase=2 skips 10,11; Codon 1 = 12..14; 15..16 is incomplete!
        ref_seq_list = ["N"] * 30
        ref_seq_list[11:14] = list("ATG")  # pos 12..14
        ref_seq_list[14:17] = list("AAA")  # pos 15..17
        ref_seq = "".join(ref_seq_list)

        gene_trailing = {
            "name": "GeneP2",
            "type": "CDS",
            "start": 10,
            "end": 16,
            "strand": "+",
            "phase": 2,
            "multi_part": False,
        }
        res_c1 = parse_metrics.get_aa_change(12, "A", "G", gene_trailing, ref_seq)
        self.assertEqual(res_c1["aa_change"], "M1V")
        self.assertIsNone(res_c1["status"])

        # pos 15 falls inside trailing 2-bp partial codon 15..16 -> must fail closed!
        res_trail = parse_metrics.get_aa_change(15, "A", "C", gene_trailing, ref_seq)
        self.assertEqual(res_trail["aa_change"], "AA annotation unavailable")
        self.assertEqual(res_trail["status"], "incomplete_partial_codon")

    def test_negative_strand_phase_1_and_2_biological_ground_truth(self):
        # CDS 10..16 (7 bases) on - strand with phase=1:
        # 5' coding base pos 16 is skipped!
        # Codon 1 is genomic 15,14,13; Codon 2 is genomic 12,11,10.
        # Put genomic 13..15 = "TAC" -> RC = "GTA" (V1)
        # Put genomic 10..12 = "CAT" -> RC = "ATG" (M2)
        ref_seq_list = ["N"] * 30
        ref_seq_list[9:12] = list("CAT")   # pos 10..12 -> RC ATG
        ref_seq_list[12:15] = list("TAC")  # pos 13..15 -> RC GTA
        ref_seq_list[15] = "A"             # pos 16 skipped (phase=1)
        ref_seq = "".join(ref_seq_list)

        gene_p1 = {
            "name": "GeneN1",
            "type": "CDS",
            "start": 10,
            "end": 16,
            "strand": "-",
            "phase": 1,
            "multi_part": False,
        }

        # pos 16 is skipped base -> fail closed
        skip_res = parse_metrics.get_aa_change(16, "A", "G", gene_p1, ref_seq)
        self.assertEqual(skip_res["aa_change"], "AA annotation unavailable")
        self.assertEqual(skip_res["status"], "incomplete_partial_codon")

        # pos 14 is A (comp T), base 2 of Codon 1 GTA -> if genomic A->C, comp G -> GGA (G1): V1G
        c1_res = parse_metrics.get_aa_change(14, "A", "C", gene_p1, ref_seq)
        self.assertEqual(c1_res["aa_change"], "V1G")
        self.assertIsNone(c1_res["status"])

        # pos 11 is A (comp T), base 2 of Codon 2 ATG -> if genomic A->G, comp C -> ACG (T2): M2T
        c2_res = parse_metrics.get_aa_change(11, "A", "G", gene_p1, ref_seq)
        self.assertEqual(c2_res["aa_change"], "M2T")
        self.assertIsNone(c2_res["status"])

        # Phase=2 on 10..17: skips 17,16; Codon 1 is 15..13 ("GTA"), pos 12..10 ("ATG") is Codon 2
        ref_seq_list[16] = "C"  # pos 17 skipped
        ref_seq2 = "".join(ref_seq_list)
        gene_p2 = {
            "name": "GeneN2",
            "type": "CDS",
            "start": 10,
            "end": 17,
            "strand": "-",
            "phase": 2,
            "multi_part": False,
        }
        skip2 = parse_metrics.get_aa_change(16, "A", "G", gene_p2, ref_seq2)
        self.assertEqual(skip2["status"], "incomplete_partial_codon")
        ok2 = parse_metrics.get_aa_change(14, "A", "C", gene_p2, ref_seq2)
        self.assertEqual(ok2["aa_change"], "V1G")

    def test_real_sars_cov2_orf1a_and_orf1ab_slippery_site(self):
        fasta_path = os.path.join(REPO_ROOT, "example_data", "illumina_SARS2", "refs", "ref.SARS-CoV-2.fasta")
        gff_path = os.path.join(REPO_ROOT, "example_data", "illumina_SARS2", "refs", "ref.SARS-CoV-2.gff")
        models_path = os.path.join(REPO_ROOT, "resources", "annotation_models", "cds_models.json")

        ref_info = parse_metrics.parse_reference(fasta_path)
        self.assertEqual(ref_info[0]["name"], "NC_045512.2")
        ref_seq = ref_info[0]["sequence"]
        genes = parse_metrics.parse_gff(gff_path)
        models = parse_metrics.load_curated_models_cached(models_path)

        cds_list = genes["NC_045512.2"]
        orf1ab_feat = next(g for g in cds_list if g.get("protein_id") == "YP_009724389.1")
        orf1a_feat = next(g for g in cds_list if g.get("protein_id") == "YP_009725295.1")

        # True SARS-CoV-2 NC_045512.2 slippery site 13462..13468 = TTTAAAC
        # Pos 13468 is 'C' (3rd base of ORF1ab codon 4401 AAC->N, AND 1st base of ORF1b codon 4402 CGG->R).
        # Mutating 13468 C->T gives: codon 4401 AAT (N4401N syn) AND codon 4402 TGG (R4402W).
        res_slip = parse_metrics.get_aa_change(13468, "C", "T", orf1ab_feat, ref_seq, curated_models=models)
        self.assertEqual(res_slip["aa_change"], "N4401N (syn)/R4402W")
        self.assertIsNone(res_slip["status"])

        # Downstream ORF1b codon 4403 is GTT (V) at 13471..13473 -> 13471 G->A => ATT (I)
        res_1b = parse_metrics.get_aa_change(13471, "G", "A", orf1ab_feat, ref_seq, curated_models=models)
        self.assertEqual(res_1b["aa_change"], "V4403I")
        self.assertEqual(res_1b["effect"], "missense")
        self.assertIsNone(res_1b["status"])

        # Contrast with colinear ORF1a at pos 13469..13471 (codon 4402 GGG -> G): 13471 G->A => GGA (G4402G syn)
        res_1a = parse_metrics.get_aa_change(13471, "G", "A", orf1a_feat, ref_seq, curated_models=models)
        self.assertEqual(res_1a["aa_change"], "G4402G (syn)")
        self.assertIsNone(res_1a["status"])

        # Reference alias support (MN908947.3) with matching sequence succeeds
        orf1ab_alias = dict(orf1ab_feat, ref_id="MN908947.3")
        res_alias = parse_metrics.get_aa_change(13471, "G", "A", orf1ab_alias, ref_seq, curated_models=models)
        self.assertEqual(res_alias["aa_change"], "V4403I")

        # Motif mismatch fails closed
        mutated_seq = ref_seq[:13461] + "GGGGGGG" + ref_seq[13468:]
        res_bad_motif = parse_metrics.get_aa_change(13471, "G", "A", orf1ab_feat, mutated_seq, curated_models=models)
        self.assertEqual(res_bad_motif["aa_change"], "AA annotation unavailable")
        self.assertEqual(res_bad_motif["status"], "curated_model_reference_mismatch")

        # 434-bp partial contig matching pattern fails closed on bounds/motif
        orf1ab_short = dict(orf1ab_feat, ref_id="SARS-CoV-2_NC_045512")
        res_short = parse_metrics.get_aa_change(100, "A", "G", orf1ab_short, "A" * 434, curated_models=models)
        self.assertEqual(res_short["aa_change"], "AA annotation unavailable")
        self.assertEqual(res_short["status"], "curated_model_reference_mismatch")

    def test_negative_strand_adjacent_minus1_frameshift(self):
        # Negative strand CDS in descending coding order:
        # Interval 1: 30..22 (9 bases: pos 30..22); Interval 2: 21..13 with frame_offset=-1 (adjacent: end 21 == prev_start 22 - 1)
        # Re-read pos 22 at start of Interval 2!
        # Let's craft genomic bases so RC coding sequence is:
        # Interval 1 RC (30->22): ATG AAA CCC (9 bases; pos 22 comp = C, genomic 22 = G)
        # Re-read pos 22 (comp C) + Interval 2 RC (21->13): C + TT TTT AAA ... => codon 4 is CTT (L)
        # Pos 22 genomic G -> A:
        #   - In Interval 1 codon 3 (pos 24,23,22 -> RC C,C,C = CCC P3), comp becomes T -> CCT (P3 syn)
        #   - In Interval 2 codon 4 (pos 22,21,20 -> RC C,T,T = CTT L4), comp becomes T -> TTT (F4)
        ref_seq_list = ["N"] * 40
        # Genomic 22..30 is RC of ATGAAACCC => revcomp("ATGAAACCC") = GGG TTT CAT
        ref_seq_list[21:30] = list("GGGTTTCAT")
        # Wait double check: idx 21 = pos 22 = 'G' (comp C); idx 29 = pos 30 = 'T' (comp A). Yes!
        # Interval 2 (21..13) RC after prepending C should have TT at 21,20 -> genomic 20..21 = AA!
        ref_seq_list[12:21] = list("TTTTTTTTA")
        ref_seq_list[19:21] = list("AA")  # pos 20..21 = A -> RC T,T
        ref_seq = "".join(ref_seq_list)

        model = {
            "model_id": "neg_adj_minus1",
            "reference": "ref_neg_fs",
            "feature_names": ["neg_fs"],
            "strand": "-",
            "translation_exception_status": "unsupported_programmed_frameshift",
            "intervals": [
                {"start": 22, "end": 30, "phase": 0, "frame_offset": 0},
                {"start": 13, "end": 21, "phase": 0, "frame_offset": -1},
            ],
        }
        parse_metrics.CURATED_CDS_MODELS_CACHE = [model]
        try:
            gene = {
                "ref_id": "ref_neg_fs",
                "name": "neg_fs",
                "strand": "-",
                "translation_exception_status": "unsupported_programmed_frameshift",
            }
            res = parse_metrics.get_aa_change(22, "G", "A", gene, ref_seq)
            self.assertEqual(res["aa_change"], "P3P (syn)/L4F")
            self.assertIsNone(res["status"])
        finally:
            parse_metrics.CURATED_CDS_MODELS_CACHE = None

    def test_plus1_frameshift_handling_and_invalid_offset_fail_closed(self):
        # Interval 1: 10..15 (ATG AAA), Interval 2: 16..21 with frame_offset=+1 skips pos 16 ('C') -> 17..19 is TTT (F3)
        ref_seq_list = ["N"] * 30
        ref_seq_list[9:15] = list("ATGAAA")
        ref_seq_list[15:21] = list("CTTTGG")  # pos 16=C skipped, 17..19=TTT (F3)
        ref_seq = "".join(ref_seq_list)

        model = {
            "model_id": "pos_plus1",
            "reference": "ref_p1",
            "feature_names": ["gene_p1"],
            "strand": "+",
            "translation_exception_status": "unsupported_programmed_frameshift",
            "intervals": [
                {"start": 10, "end": 15, "phase": 0, "frame_offset": 0},
                {"start": 16, "end": 21, "phase": 0, "frame_offset": 1},
            ],
        }
        parse_metrics.CURATED_CDS_MODELS_CACHE = [model]
        try:
            gene = {
                "ref_id": "ref_p1",
                "name": "gene_p1",
                "strand": "+",
                "translation_exception_status": "unsupported_programmed_frameshift",
            }
            # Skipped base 16 is out of translated CDS bounds
            res_skip = parse_metrics.get_aa_change(16, "C", "A", gene, ref_seq)
            self.assertEqual(res_skip["aa_change"], "AA annotation unavailable")
            self.assertEqual(res_skip["status"], "curated_model_position_out_of_bounds")

            # Pos 17 is first base of codon 3 TTT -> ATT (F3I)
            res_c3 = parse_metrics.get_aa_change(17, "T", "A", gene, ref_seq)
            self.assertEqual(res_c3["aa_change"], "F3I")
            self.assertIsNone(res_c3["status"])

            # Invalid frame_offset=2 fails closed
            model["intervals"][1]["frame_offset"] = 2
            res_bad = parse_metrics.get_aa_change(17, "T", "A", gene, ref_seq)
            self.assertEqual(res_bad["aa_change"], "AA annotation unavailable")
            self.assertEqual(res_bad["status"], "curated_model_boundary_unsupported")
        finally:
            parse_metrics.CURATED_CDS_MODELS_CACHE = None

    def test_coding_indel_classification(self):
        gene = {"name": "S", "type": "CDS", "start": 100, "end": 500, "strand": "+", "phase": 0, "multi_part": False}
        ref_seq = "A" * 600
        del_9 = parse_metrics.get_aa_change(120, "GTCTGGTTTT", "G", gene, ref_seq)
        self.assertEqual(del_9["aa_change"], "In-frame del (3 aa)")
        self.assertEqual(del_9["effect"], "inframe_deletion")

        ins_3 = parse_metrics.get_aa_change(120, "A", "ACGT"[:4], gene, ref_seq)
        self.assertEqual(ins_3["aa_change"], "In-frame ins (1 aa)")
        self.assertEqual(ins_3["effect"], "inframe_insertion")

        fs_ins = parse_metrics.get_aa_change(120, "C", "CT", gene, ref_seq)
        self.assertEqual(fs_ins["aa_change"], "Frameshift ins (+1 bp)")
        self.assertEqual(fs_ins["effect"], "frameshift_insertion")

        fs_del = parse_metrics.get_aa_change(120, "ATG", "A", gene, ref_seq)
        self.assertEqual(fs_del["aa_change"], "Frameshift del (-2 bp)")
        self.assertEqual(fs_del["effect"], "frameshift_deletion")

    def test_multi_part_cds_detection_and_bypass(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "multipart.gff")
            with open(gff_path, "w") as f:
                f.writelines([
                    "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tParent=gene_1\n",
                    "ref_contig\tSource\tCDS\t200\t250\t.\t+\t1\tParent=gene_1\n",
                ])
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 2)
            for g in ref_genes:
                self.assertTrue(g["multi_part"])
                res = parse_metrics.get_aa_change(120, "A", "G", g, "A" * 300)
                self.assertEqual(res["aa_change"], "AA annotation unavailable")
                self.assertEqual(res["status"], "unsupported_multipart_cds")

    def test_weak_names_and_unknowns_not_grouped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "weak.gff")
            with open(gff_path, "w") as f:
                f.writelines([
                    "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tgene=WeakSharedName\n",
                    "ref_contig\tSource\tCDS\t200\t250\t.\t+\t0\tgene=WeakSharedName\n",
                    "ref_contig\tSource\tCDS\t300\t350\t.\t+\t0\tParent=Unknown\n",
                    "ref_contig\tSource\tCDS\t400\t450\t.\t+\t0\tParent=Unknown\n",
                ])
            genes = parse_metrics.parse_gff(gff_path)
            for g in genes["ref_contig"]:
                self.assertFalse(g["multi_part"])

    def test_translation_exception_flags_and_underscore_normalization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "exc.gff")
            with open(gff_path, "w") as f:
                f.writelines([
                    "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tID=f1;exception=ribosomal_slippage\n",
                    "ref_contig\tSource\tCDS\t200\t250\t.\t+\t0\tID=f2;Note=programmed_ribosomal_frameshift\n",
                    "ref_contig\tSource\tCDS\t300\t350\t.\t+\t0\tID=f3;transl_except=(pos:310..312%2Caa:Sec)\n",
                ])
            genes = parse_metrics.parse_gff(gff_path)["ref_contig"]
            self.assertEqual(genes[0]["translation_exception_status"], "unsupported_programmed_frameshift")
            self.assertTrue(genes[0]["requires_curated_model"])
            self.assertEqual(genes[1]["translation_exception_status"], "unsupported_programmed_frameshift")
            self.assertEqual(genes[2]["translation_exception_status"], "unsupported_translation_exception")

    def test_cds_biological_sorting_and_distinct_protein_grouping(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "sort.gff")
            with open(gff_path, "w") as f:
                f.writelines([
                    "ref_contig\tSource\tCDS\t200\t250\t.\t+\t0\tParent=gene_pos\n",
                    "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tParent=gene_pos\n",
                    "ref_contig\tSource\tCDS\t300\t350\t.\t-\t0\tParent=gene_neg\n",
                    "ref_contig\tSource\tCDS\t400\t450\t.\t-\t0\tParent=gene_neg\n",
                    "ref_contig\tRefSeq\tCDS\t500\t600\t.\t+\t0\tID=cds-ab;Parent=gene-1;protein_id=YP_ab;exception=ribosomal slippage\n",
                    "ref_contig\tRefSeq\tCDS\t600\t800\t.\t+\t0\tID=cds-ab;Parent=gene-1;protein_id=YP_ab;exception=ribosomal slippage\n",
                    "ref_contig\tRefSeq\tCDS\t500\t620\t.\t+\t0\tID=cds-a;Parent=gene-1;protein_id=YP_a;product=ORF1a polyprotein\n",
                ])
            genes = parse_metrics.parse_gff(gff_path)["ref_contig"]
            pos_parts = [g for g in genes if g["parent"] == "gene_pos"]
            neg_parts = [g for g in genes if g["parent"] == "gene_neg"]
            self.assertEqual([p["start"] for p in pos_parts], [100, 200])
            self.assertEqual([p["start"] for p in neg_parts], [400, 300])
            orf1a = [g for g in genes if g["protein_id"] == "YP_a"]
            orf1ab = [g for g in genes if g["protein_id"] == "YP_ab"]
            self.assertFalse(orf1a[0]["multi_part"])
            self.assertEqual(len(orf1ab), 2)
            self.assertTrue(all(g["multi_part"] for g in orf1ab))

    def test_overlapping_genes_structured_consequences_no_garbled_join(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            gff_path = os.path.join(tmpdir, "overlap.gff")
            with open(gff_path, "w") as f:
                f.writelines([
                    "ref1\tSrc\tCDS\t100\t200\t.\t+\t0\tID=ga;gene=GeneA;protein_id=P_A;product=Protein A\n",
                    "ref1\tSrc\tCDS\t180\t280\t.\t+\t0\tID=gb;gene=GeneB;protein_id=P_B;exception=ribosomal slippage\n",
                ])
            vcf_path = os.path.join(tmpdir, "overlap.vcf")
            with open(vcf_path, "w") as f:
                f.writelines([
                    "##fileformat=VCFv4.2\n",
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n",
                    "ref1\t190\t.\tA\tC\t50\tPASS\tVAF=1.0;CONSENSUS\tDP:AD\t100:0,100\n",
                ])
            ref_seq = "A" * 300
            genes = parse_metrics.parse_gff(gff_path)
            vars_res = parse_metrics.parse_vcf(vcf_path, genes, ref_seqs={"ref1": ref_seq}, curated_models=[])
            self.assertEqual(len(vars_res), 1)
            v = vars_res[0]
            self.assertEqual(v["gene"], "GeneA/GeneB")
            # Primary summary aa_change should NOT append '/AA annotation unavailable' to resolved GeneA annotation!
            self.assertNotIn("AA annotation unavailable", v["aa_change"])
            self.assertEqual(len(v["consequences"]), 2)
            self.assertEqual(v["consequences"][0]["gene"], "GeneA")
            self.assertIsNone(v["consequences"][0]["status"])
            self.assertEqual(v["consequences"][1]["gene"], "GeneB")
            self.assertEqual(v["consequences"][1]["aa_change"], "AA annotation unavailable")
            self.assertEqual(v["consequences"][1]["status"], "unsupported_programmed_frameshift")


if __name__ == "__main__":
    unittest.main()
