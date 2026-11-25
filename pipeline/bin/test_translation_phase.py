#!/usr/bin/env python3
import unittest
import sys
import os

# Adjust path to import parse_metrics
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import parse_metrics

class TestTranslationPhase(unittest.TestCase):
    
    def test_positive_strand_phase_1(self):
        ref_seq_list = ["N"] * 200
        ref_seq_list[98] = "A"
        ref_seq_list[99] = "T"
        ref_seq_list[100] = "G"
        ref_seq_list[101] = "A"
        ref_seq_list[102] = "C" # mutation at coordinate 103 (index 102) C -> G (alt)
        ref_seq_list[103] = "C"
        ref_seq = "".join(ref_seq_list)
        
        gene = {
            "name": "GeneP",
            "type": "CDS",
            "start": 100,
            "end": 150,
            "strand": "+",
            "phase": 1,
            "multi_part": False
        }
        
        res = parse_metrics.get_aa_change(103, "C", "G", gene, ref_seq)
        self.assertEqual(res["aa_change"], "T2S", f"Expected positive-strand phase=1 mutation to yield T2S, got: {res}")
        self.assertIsNone(res["status"])
        
    def test_negative_strand_phase_1(self):
        ref_seq_list = ["N"] * 50
        ref_seq_list[10] = "G"
        ref_seq_list[11] = "T"
        ref_seq_list[12] = "A"
        ref_seq_list[13] = "C"
        ref_seq_list[14] = "G"
        ref_seq_list[15] = "T"
        ref_seq = "".join(ref_seq_list)
        
        gene = {
            "name": "GeneN",
            "type": "CDS",
            "start": 10,
            "end": 15,
            "strand": "-",
            "phase": 1,
            "multi_part": False
        }
        
        res = parse_metrics.get_aa_change(12, "T", "C", gene, ref_seq)
        self.assertEqual(res["aa_change"], "Y2C", f"Expected negative-strand phase=1 mutation to yield Y2C, got: {res}")
        self.assertIsNone(res["status"])

    def test_multi_part_cds_detection_and_bypass(self):
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tParent=gene_1\n",
            "ref_contig\tSource\tCDS\t200\t250\t.\t+\t1\tParent=gene_1\n"
        ]
        
        gff_path = "temp_test_multipart.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
            
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            
            self.assertEqual(len(ref_genes), 2)
            for g in ref_genes:
                self.assertTrue(g["multi_part"], "Spliced CDS sharing strong Parent should be tagged as multi_part")
                
                res = parse_metrics.get_aa_change(120, "A", "G", g, "N"*300)
                self.assertEqual(res["aa_change"], "AA annotation unavailable", "Should bypass translation for multi-part CDS")
                self.assertIsNone(res["status"])
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_weak_names_not_grouped(self):
        # Unrelated features with the same weak name/gene are not falsely marked multi-part
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tgene=WeakSharedName\n",
            "ref_contig\tSource\tCDS\t200\t250\t.\t+\t0\tgene=WeakSharedName\n",
            "ref_contig\tSource\tCDS\t300\t350\t.\t+\t0\tName=WeakSharedName\n",
            "ref_contig\tSource\tCDS\t400\t450\t.\t+\t0\tproduct=WeakSharedName\n"
        ]
        gff_path = "temp_test_weaknames.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 4)
            for g in ref_genes:
                self.assertFalse(g["multi_part"], f"Features grouped by weak name alone! {g}")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_unknown_features_not_grouped(self):
        # Unknown features are not grouped together solely because their name/value is Unknown
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tParent=Unknown\n",
            "ref_contig\tSource\tCDS\t200\t250\t.\t+\t0\tParent=Unknown\n",
            "ref_contig\tSource\tCDS\t300\t350\t.\t+\t0\tID=Unknown\n",
            "ref_contig\tSource\tCDS\t400\t450\t.\t+\t0\tgene=Unknown\n"
        ]
        gff_path = "temp_test_unknowns.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 4)
            for g in ref_genes:
                self.assertFalse(g["multi_part"], f"Unknown features falsely grouped together! {g}")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_programmed_frameshift_exception(self):
        # a CDS with exception=programmed ribosomal frameshift returns AA annotation unavailable
        # with status unsupported_programmed_frameshift
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tID=ORF1ab;exception=programmed ribosomal frameshift\n"
        ]
        gff_path = "temp_test_frameshift.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 1)
            gene = ref_genes[0]
            self.assertEqual(gene["translation_exception_status"], "unsupported_programmed_frameshift")
            
            res = parse_metrics.get_aa_change(120, "A", "G", gene, "N"*300)
            self.assertEqual(res["aa_change"], "AA annotation unavailable")
            self.assertEqual(res["status"], "unsupported_programmed_frameshift")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_ribosomal_slippage_note(self):
        # a CDS with Note mentioning ribosomal slippage returns AA annotation unavailable
        # with status unsupported_programmed_frameshift
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tID=ORF1ab;Note=ribosomal slippage occurs here\n"
        ]
        gff_path = "temp_test_slippage.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 1)
            gene = ref_genes[0]
            self.assertEqual(gene["translation_exception_status"], "unsupported_programmed_frameshift")
            
            res = parse_metrics.get_aa_change(120, "A", "G", gene, "N"*300)
            self.assertEqual(res["aa_change"], "AA annotation unavailable")
            self.assertEqual(res["status"], "unsupported_programmed_frameshift")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_transl_except_status(self):
        # transl_except is present without explicit frameshift/slippage wording,
        # set translation_exception_status: "unsupported_translation_exception"
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tID=gene_te;transl_except=(pos:110..112%2Caa:Sec)\n"
        ]
        gff_path = "temp_test_te.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 1)
            gene = ref_genes[0]
            self.assertEqual(gene["translation_exception_status"], "unsupported_translation_exception")
            
            res = parse_metrics.get_aa_change(120, "A", "G", gene, "N"*300)
            self.assertEqual(res["aa_change"], "AA annotation unavailable")
            self.assertEqual(res["status"], "unsupported_translation_exception")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_underscore_normalization_slippage(self):
        # Test exception=ribosomal_slippage is normalized and sets requires_curated_model = True
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tID=ORF1ab;exception=ribosomal_slippage\n"
        ]
        gff_path = "temp_test_norm_slippage.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 1)
            gene = ref_genes[0]
            self.assertEqual(gene["translation_exception_status"], "unsupported_programmed_frameshift")
            self.assertTrue(gene["requires_curated_model"])
            
            res = parse_metrics.get_aa_change(120, "A", "G", gene, "N"*300)
            self.assertEqual(res["aa_change"], "AA annotation unavailable")
            self.assertEqual(res["status"], "unsupported_programmed_frameshift")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_underscore_normalization_frameshift(self):
        # Test exception=programmed_ribosomal_frameshift is normalized and sets requires_curated_model = True
        gff_lines = [
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tID=ORF1ab;exception=programmed_ribosomal_frameshift\n"
        ]
        gff_path = "temp_test_norm_frameshift.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 1)
            gene = ref_genes[0]
            self.assertEqual(gene["translation_exception_status"], "unsupported_programmed_frameshift")
            self.assertTrue(gene["requires_curated_model"])
            
            res = parse_metrics.get_aa_change(120, "A", "G", gene, "N"*300)
            self.assertEqual(res["aa_change"], "AA annotation unavailable")
            self.assertEqual(res["status"], "unsupported_programmed_frameshift")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_cds_biological_sorting(self):
        # Sort CDS parts in biological coding order:
        # + strand: ascending start, then ascending end
        # - strand: descending start, then descending end
        # We will check the order of sorted CDS parts returned by parse_gff
        gff_lines = [
            # + strand out of order in GFF
            "ref_contig\tSource\tCDS\t200\t250\t.\t+\t0\tParent=gene_pos\n",
            "ref_contig\tSource\tCDS\t100\t150\t.\t+\t0\tParent=gene_pos\n",
            # - strand out of order in GFF
            "ref_contig\tSource\tCDS\t300\t350\t.\t-\t0\tParent=gene_neg\n",
            "ref_contig\tSource\tCDS\t400\t450\t.\t-\t0\tParent=gene_neg\n"
        ]
        gff_path = "temp_test_sort.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 4)
            
            # Since global sort is by start position ascending, let's look at the grouping keys
            pos_parts = [g for g in ref_genes if g["parent"] == "gene_pos"]
            neg_parts = [g for g in ref_genes if g["parent"] == "gene_neg"]
            
            # Verify + strand parts order (ascending start)
            self.assertEqual(pos_parts[0]["start"], 100)
            self.assertEqual(pos_parts[1]["start"], 200)
            
            # Verify - strand parts order (descending start)
            self.assertEqual(neg_parts[0]["start"], 400)
            self.assertEqual(neg_parts[1]["start"], 300)
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_toy_coronavirus_curated_model(self):
        # Setup toy reference sequence
        ref_seq_list = ["N"] * 50
        # First interval: 10..19. Let sequence be ATG AAA CCC G (length 10)
        # Positions: 10, 11, 12, 13, 14, 15, 16, 17, 18, 19
        # Second interval: 19..30. Let sequence be G AC TTT AAA GG (length 12)
        # Overlap is at 19 (which is G).
        # Stitched CDS: ATG AAA CCC G AC TTT AAA GG (length 22, 7 codons + 1 base)
        # Sequence at indices:
        # 10(A), 11(T), 12(G) -> ATG (M)
        # 13(A), 14(A), 15(A) -> AAA (K)
        # 16(C), 17(C), 18(C) -> CCC (P)
        # 19(G), 20(A), 21(C) -> GAC (D)
        # 22(T), 23(T), 24(T) -> TTT (F)
        # ...
        ref_seq_list[9] = "A"
        ref_seq_list[10] = "T"
        ref_seq_list[11] = "G"
        ref_seq_list[12] = "A"
        ref_seq_list[13] = "A"
        ref_seq_list[14] = "A"
        ref_seq_list[15] = "C"
        ref_seq_list[16] = "C"
        ref_seq_list[17] = "C"
        ref_seq_list[18] = "G"
        ref_seq_list[19] = "A"
        ref_seq_list[20] = "C"
        ref_seq_list[21] = "T"
        ref_seq_list[22] = "T"
        ref_seq_list[23] = "T"
        ref_seq = "".join(ref_seq_list)

        toy_model = {
            "model_id": "toy_pos",
            "reference": "toy_ref",
            "feature_names": ["toy_gene"],
            "strand": "+",
            "translation_exception_status": "unsupported_programmed_frameshift",
            "intervals": [
                {"start": 10, "end": 19, "phase": 0},
                {"start": 19, "end": 30, "phase": 0}
            ]
        }
        
        # Override cache
        parse_metrics.CURATED_CDS_MODELS_CACHE = [toy_model]
        try:
            gene = {
                "ref_id": "toy_ref",
                "name": "toy_gene",
                "strand": "+",
                "translation_exception_status": "unsupported_programmed_frameshift"
            }
            # Mutation at pos 14: A -> C (second base of codon 2 AAA -> ACA, translates to T)
            res = parse_metrics.get_aa_change(14, "A", "C", gene, ref_seq)
            self.assertEqual(res["aa_change"], "K2T")
            self.assertIsNone(res["status"])
            
            # Mutation at pos 19 (overlap, represented twice): G -> A
            # Under overlap, maps to index 9 (first base of codon 4 GGA -> AGA, R) and index 10 (second base of codon 4 GGA -> GAA, E)
            res_overlap = parse_metrics.get_aa_change(19, "G", "A", gene, ref_seq)
            self.assertEqual(res_overlap["aa_change"], "G4R/G4E")
        finally:
            parse_metrics.CURATED_CDS_MODELS_CACHE = None

    def test_toy_coronavirus_curated_model_negative_strand(self):
        # Test generic negative-strand curated model handling
        # Reference sequence in reverse direction:
        # CDS (genomic 15..20): 20(T), 19(C), 18(A), 17(C), 16(A), 15(T)
        # Reverse complement:
        # 15(T->A), 16(A->T), 17(C->G), 18(A->T), 19(C->G), 20(T->A)
        # Combined RC sequence: AGT GTA (translates to S, V)
        ref_seq_list = ["N"] * 30
        ref_seq_list[14] = "T" # pos 15
        ref_seq_list[15] = "A" # pos 16
        ref_seq_list[16] = "C" # pos 17
        ref_seq_list[17] = "A" # pos 18
        ref_seq_list[18] = "C" # pos 19
        ref_seq_list[19] = "T" # pos 20
        ref_seq = "".join(ref_seq_list)

        toy_model = {
            "model_id": "toy_neg",
            "reference": "toy_ref_neg",
            "feature_names": ["toy_gene_neg"],
            "strand": "-",
            "translation_exception_status": "unsupported_programmed_frameshift",
            "intervals": [
                {"start": 15, "end": 20, "phase": 0}
            ]
        }
        
        parse_metrics.CURATED_CDS_MODELS_CACHE = [toy_model]
        try:
            gene = {
                "ref_id": "toy_ref_neg",
                "name": "toy_gene_neg",
                "strand": "-",
                "translation_exception_status": "unsupported_programmed_frameshift"
            }
            # Mutation at pos 17 (genomic C). In RC sequence, codon 2 is GTA.
            # Genomic pos 17 is complement G (first base of codon 2). Mutate genomic C -> A (RC is T).
            # Codon becomes TTA (L).
            # Expected change: V2L
            res = parse_metrics.get_aa_change(17, "C", "A", gene, ref_seq)
            self.assertEqual(res["aa_change"], "V2L")
        finally:
            parse_metrics.CURATED_CDS_MODELS_CACHE = None

    def test_unmatched_frameshift_fails_closed(self):
        # negative test showing unmatched frameshift/slippage features still fail closed
        gene = {
            "ref_id": "unknown_ref",
            "name": "unknown_gene",
            "strand": "+",
            "translation_exception_status": "unsupported_programmed_frameshift"
        }
        parse_metrics.CURATED_CDS_MODELS_CACHE = []
        try:
            res = parse_metrics.get_aa_change(120, "A", "G", gene, "N"*300)
            self.assertEqual(res["aa_change"], "AA annotation unavailable")
            self.assertEqual(res["status"], "unsupported_programmed_frameshift")
        finally:
            parse_metrics.CURATED_CDS_MODELS_CACHE = None

    def test_real_sars_cov2_orf1ab_coordinate(self):
        # Test real or miniature NC_045512.2-style ORF1ab coordinate test with a hand-checked expected AA change
        # Slippery site: 13462..13468
        # Sequence at 13462-13468 is TTT AAA C
        # Codon 4401: index 13466(A), 13467(A), 13468(A) -> AAA (K)
        # Let's verify!
        # First interval ends at 13468. Second interval starts at 13468.
        # Let's load the real model from the config file itself to ensure everything matches:
        models = parse_metrics.load_curated_models_cached()
        self.assertTrue(len(models) > 0, "Curated models file should not be empty")
        
        # Mock reference sequence spanning the target area
        # We need a sequence of at least 21555 bases for the NC_045512.2 reference
        ref_seq_list = ["N"] * 22000
        # Set codon 4401 at coordinates 13466, 13467, 13468 (indices 13465, 13466, 13467) to A, A, A
        ref_seq_list[13465] = "A"
        ref_seq_list[13466] = "A"
        ref_seq_list[13467] = "A"
        ref_seq = "".join(ref_seq_list)
        
        gene = {
            "ref_id": "NC_045512.2",
            "name": "ORF1ab",
            "strand": "+",
            "translation_exception_status": "unsupported_programmed_frameshift"
        }
        
        # Mutation at coordinate 13467 (index 13466) A -> G.
        # This is offset 1 of codon 4401 (AAA -> AGA, translates to R)
        # Expected change: K4401R
        res = parse_metrics.get_aa_change(13467, "A", "G", gene, ref_seq)
        self.assertIn("K4401R", res["aa_change"])
        self.assertIsNone(res["status"])

    def test_distinct_proteins_sharing_parent_not_collapsed(self):
        # RefSeq-style annotation where ORF1a (pp1a) and ORF1ab (pp1ab) share Parent=gene-1,
        # but have distinct protein_ids (YP_009724389.1 vs YP_009725295.1).
        # ORF1ab has ribosomal slippage; ORF1a is a colinear single-part CDS.
        gff_lines = [
            "ref_contig\tRefSeq\tCDS\t266\t13468\t.\t+\t0\tID=cds-ab;Parent=gene-1;protein_id=YP_ab;exception=ribosomal slippage\n",
            "ref_contig\tRefSeq\tCDS\t13468\t21555\t.\t+\t0\tID=cds-ab;Parent=gene-1;protein_id=YP_ab;exception=ribosomal slippage\n",
            "ref_contig\tRefSeq\tCDS\t266\t13483\t.\t+\t0\tID=cds-a;Parent=gene-1;protein_id=YP_a;product=ORF1a polyprotein\n"
        ]
        gff_path = "temp_test_distinct_parent.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
        try:
            genes = parse_metrics.parse_gff(gff_path)
            ref_genes = genes.get("ref_contig", [])
            self.assertEqual(len(ref_genes), 3)
            
            orf1a_feats = [g for g in ref_genes if g["protein_id"] == "YP_a"]
            orf1ab_feats = [g for g in ref_genes if g["protein_id"] == "YP_ab"]
            
            self.assertEqual(len(orf1a_feats), 1)
            self.assertFalse(orf1a_feats[0]["multi_part"])
            self.assertIsNone(orf1a_feats[0]["translation_exception_status"])
            
            self.assertEqual(len(orf1ab_feats), 2)
            for g in orf1ab_feats:
                self.assertTrue(g["multi_part"])
                self.assertEqual(g["translation_exception_status"], "unsupported_programmed_frameshift")
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)

    def test_overlapping_reading_frames_in_parse_vcf(self):
        # Create GFF with two overlapping genes: GeneA (100..200) and GeneB (180..280)
        gff_lines = [
            "ref1\tSrc\tCDS\t100\t200\t.\t+\t0\tID=ga;gene=GeneA\n",
            "ref1\tSrc\tCDS\t180\t280\t.\t+\t0\tID=gb;gene=GeneB\n"
        ]
        gff_path = "temp_test_overlap.gff"
        with open(gff_path, "w") as f:
            f.writelines(gff_lines)
            
        # Create dummy VCF with a variant at 190 (which is in both GeneA and GeneB)
        vcf_lines = [
            "##fileformat=VCFv4.2\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            "ref1\t190\t.\tA\tC\t50\tPASS\tDP=100;AF=1.0\n"
        ]
        vcf_path = "temp_test_overlap.vcf"
        with open(vcf_path, "w") as f:
            f.writelines(vcf_lines)
            
        ref_seq = "A" * 300
        try:
            genes = parse_metrics.parse_gff(gff_path)
            vars_res = parse_metrics.parse_vcf(vcf_path, genes, ref_seqs={"ref1": ref_seq})
            self.assertEqual(len(vars_res), 1)
            v = vars_res[0]
            self.assertEqual(v["gene"], "GeneA/GeneB")
            self.assertIn("/", v["aa_change"])
        finally:
            if os.path.exists(gff_path):
                os.remove(gff_path)
            if os.path.exists(vcf_path):
                os.remove(vcf_path)

    def test_frame_offset_adjacent_intervals(self):
        ref_seq_list = ["N"] * 50
        ref_seq_list[9] = "A"; ref_seq_list[10] = "T"; ref_seq_list[11] = "G"
        ref_seq_list[12] = "A"; ref_seq_list[13] = "A"; ref_seq_list[14] = "A"
        ref_seq_list[15] = "C"; ref_seq_list[16] = "C"; ref_seq_list[17] = "C" # 18: C
        ref_seq_list[18] = "A"; ref_seq_list[19] = "A"                         # 19, 20: A
        ref_seq = "".join(ref_seq_list)

        model = {
            "model_id": "adj_model",
            "reference": "ref_adj",
            "feature_names": ["adj_gene"],
            "strand": "+",
            "translation_exception_status": "unsupported_programmed_frameshift",
            "intervals": [
                {"start": 10, "end": 18, "phase": 0, "frame_offset": 0},
                {"start": 19, "end": 28, "phase": 0, "frame_offset": -1}
            ]
        }
        parse_metrics.CURATED_CDS_MODELS_CACHE = [model]
        try:
            gene = {
                "ref_id": "ref_adj",
                "name": "adj_gene",
                "strand": "+",
                "translation_exception_status": "unsupported_programmed_frameshift"
            }
            res = parse_metrics.get_aa_change(18, "C", "T", gene, ref_seq)
            self.assertIn("/", res["aa_change"])
            self.assertIsNone(res["status"])
        finally:
            parse_metrics.CURATED_CDS_MODELS_CACHE = None

if __name__ == "__main__":
    unittest.main()
