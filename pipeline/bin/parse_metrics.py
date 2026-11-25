#!/usr/bin/env python3
import os
import sys
import json
import argparse
import subprocess
import urllib.parse
from collections import defaultdict

CURATED_CDS_MODELS_CACHE = None

def load_curated_models_cached():
    global CURATED_CDS_MODELS_CACHE
    if CURATED_CDS_MODELS_CACHE is not None:
        return CURATED_CDS_MODELS_CACHE
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    paths = [
        os.path.join(script_dir, "resources", "annotation_models", "cds_models.json"),
        os.path.join(script_dir, "..", "..", "resources", "annotation_models", "cds_models.json"),
        # Also check current working directory
        "resources/annotation_models/cds_models.json"
    ]
    for p in paths:
        if os.path.exists(p):
            try:
                with open(p, "r") as f:
                    CURATED_CDS_MODELS_CACHE = json.load(f)
                    return CURATED_CDS_MODELS_CACHE
            except Exception as e:
                print(f"Warning: Failed to parse curated models file {p}: {e}")
    
    CURATED_CDS_MODELS_CACHE = []
    return CURATED_CDS_MODELS_CACHE


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate NGS run metrics for Viromics Dashboard")
    parser.add_argument("--run-id", required=True, help="Unique identifier for the run")
    parser.add_argument("--fastp", required=True, help="Comma-separated list of fastp JSON files")
    parser.add_argument("--bams", required=True, help="Comma-separated list of BAM files")
    parser.add_argument("--vcfs", required=True, help="Comma-separated list of VCF files")
    parser.add_argument("--consensuses", required=True, help="Comma-separated list of consensus FASTA files")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--gff", help="Optional GFF3 annotation file")
    parser.add_argument("--primers", help="Optional path to primer scheme BED file")
    parser.add_argument("--amplicons-bed", help="Optional path to explicit amplicon coordinates BED file")
    parser.add_argument("--metadata", help="Optional path to sample metadata CSV file")
    parser.add_argument("--min-coverage", type=int, default=10, help="Minimum depth for coverage calculations")
    parser.add_argument("--min-freq", type=float, default=0.2, help="Minimum alternative allele frequency for variants")
    parser.add_argument("--output", default="viromics_run_report.json", help="Output JSON path")
    return parser.parse_args()

def parse_reference(ref_path):
    """Parse all reference sequences from a multi-FASTA file."""
    references = []
    current_name = None
    current_seq = []
    with open(ref_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if current_name:
                    seq = "".join(current_seq)
                    references.append({
                        "name": current_name,
                        "length": len(seq),
                        "sequence": seq
                    })
                current_name = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_name:
            seq = "".join(current_seq)
            references.append({
                "name": current_name,
                "length": len(seq),
                "sequence": seq
            })
    return references

def parse_gff(gff_path):
    """Parse genes and CDS from GFF3 file, grouped by reference sequence ID."""
    features = defaultdict(list)
    if not gff_path or not os.path.exists(gff_path):
        return dict(features)
        
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
                
            ref_id = parts[0]
            feature_type = parts[2]
            if feature_type in ["gene", "CDS", "exon"]:
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                
                phase = 0
                if feature_type == "CDS" and len(parts) >= 8:
                    phase_col = parts[7].strip()
                    if phase_col in [".", "", "null"]:
                        phase = 0
                    elif phase_col in ["0", "1", "2"]:
                        phase = int(phase_col)
                    else:
                        phase = -1 # Invalid phase
                
                # Parse attributes
                attrs = {}
                for attr in parts[8].split(";"):
                    if not attr.strip():
                        continue
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        attrs[k.strip()] = urllib.parse.unquote(v.strip())
                
                name = attrs.get("gene", attrs.get("Name", attrs.get("ID", "Unknown")))
                product = attrs.get("product", "")
                
                features[ref_id].append({
                    "ref_id": ref_id,
                    "name": name,
                    "type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "phase": phase,
                    "product": product,
                    "id": attrs.get("ID", ""),
                    "parent": attrs.get("Parent", ""),
                    "gene": attrs.get("gene", ""),
                    "feature_name": attrs.get("Name", ""),
                    "locus_tag": attrs.get("locus_tag", ""),
                    "protein_id": attrs.get("protein_id", ""),
                    "transcript_id": attrs.get("transcript_id", ""),
                    "exception": attrs.get("exception", ""),
                    "note": attrs.get("Note", ""),
                    "transl_except": attrs.get("transl_except", ""),
                    "attributes": attrs
                })
    
    # Hierarchical filtering per reference: CDS > gene > exon to prevent duplication clutter
    sorted_genes = {}
    for ref_id, f_list in features.items():
        # Check if CDS exists
        cds_features = [f for f in f_list if f["type"] == "CDS"]
        if cds_features:
            selected = cds_features
        else:
            gene_features = [f for f in f_list if f["type"] == "gene"]
            if gene_features:
                selected = gene_features
            else:
                selected = f_list
                
        # Separate CDS features from non-CDS features to handle grouping and sorting
        cds_only = [f for f in selected if f["type"] == "CDS"]
        non_cds_only = [f for f in selected if f["type"] != "CDS"]
        
        # Group CDS features by (ref_id, strand, grouping_key)
        cds_groups = defaultdict(list)
        ungrouped_count = 0
        for f in cds_only:
            attrs = f["attributes"]
            grouping_key = None
            for key in ["protein_id", "transcript_id", "ID", "Parent", "locus_tag"]:
                val = attrs.get(key, "").strip()
                if val and val.lower() != "unknown":
                    grouping_key = val
                    break
            
            if not grouping_key:
                ungrouped_count += 1
                grouping_key = f"__ungrouped_cds_{ungrouped_count}__"
                
            f["grouping_key"] = grouping_key
            cds_groups[(f["strand"], grouping_key)].append(f)
            
        # Sort and flag CDS parts within each group
        groups_list = []
        for (strand, gkey), group_feats in cds_groups.items():
            # Sort CDS parts in biological coding order:
            # + strand: ascending start, then ascending end
            # - strand: descending start, then descending end
            if strand == "+":
                group_feats.sort(key=lambda x: (x["start"], x["end"]))
            else:
                group_feats.sort(key=lambda x: (-x["start"], -x["end"]))
                
            # Determine multi_part status (two or more rows with distinct intervals)
            distinct_intervals = {(f["start"], f["end"]) for f in group_feats}
            is_multi = (len(group_feats) >= 2 and len(distinct_intervals) >= 2)
            
            # Determine translation_exception_status
            has_fs = False
            has_te = False
            for f in group_feats:
                exc_val = " ".join(f.get("exception", "").lower().replace("_", " ").split())
                note_val = " ".join(f.get("note", "").lower().replace("_", " ").split())
                if ("programmed ribosomal frameshift" in exc_val or
                    "ribosomal slippage" in exc_val or
                    "frameshift" in exc_val or
                    "programmed ribosomal frameshift" in note_val or
                    "ribosomal slippage" in note_val or
                    "frameshift" in note_val):
                    has_fs = True
                if f.get("transl_except", "").strip():
                    has_te = True
                    
            if has_fs:
                exception_status = "unsupported_programmed_frameshift"
            elif has_te:
                exception_status = "unsupported_translation_exception"
            else:
                exception_status = None
                
            requires_curated = (exception_status == "unsupported_programmed_frameshift")
                
            for f in group_feats:
                f["multi_part"] = is_multi
                f["translation_exception_status"] = exception_status
                f["requires_curated_model"] = requires_curated
                
            groups_list.append(group_feats)
            
        # Set default status/multi_part for non-CDS features
        for f in non_cds_only:
            f["multi_part"] = False
            f["translation_exception_status"] = None
            f["requires_curated_model"] = False
            f["grouping_key"] = None
            groups_list.append([f])
            
        # Sort the list of groups based on the minimum start position of the group
        groups_list.sort(key=lambda g: min(f["start"] for f in g))
        
        # Flatten the sorted groups into the final list
        sorted_genes[ref_id] = [f for g in groups_list for f in g]
        
    return sorted_genes

def parse_fastp(fastp_path):
    """Parse quality stats from fastp JSON report."""
    if not os.path.exists(fastp_path):
        return {}
    try:
        with open(fastp_path, 'r') as f:
            data = json.load(f)
        
        summary = data.get("summary", {})
        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})
        
        return {
            "raw_reads": before.get("total_reads", 0),
            "raw_bases": before.get("total_bases", 0),
            "trimmed_reads": after.get("total_reads", 0),
            "trimmed_bases": after.get("total_bases", 0),
            "q30_rate_raw": before.get("q30_rate", 0.0),
            "q30_rate_trimmed": after.get("q30_rate", 0.0),
            "duplication_rate": summary.get("duplication", {}).get("rate", 0.0)
        }
    except Exception as e:
        print(f"Error parsing fastp JSON {fastp_path}: {e}")
        return {}

def calculate_coverage(bam_path, ref_lengths, num_bins=500, min_cov=10):
    """Calculate depth of coverage for each reference sequence and downsample for visualization."""
    results = {}
    
    # Initialize depths list and mapped reads counts for each reference
    depths = {ref_name: [0] * ref_len for ref_name, ref_len in ref_lengths.items()}
    mapped_reads = {ref_name: 0 for ref_name in ref_lengths}
    
    # Run samtools view -c -F 2308 to retrieve primary mapped reads per contig
    for ref_name in ref_lengths:
        try:
            cmd = ["samtools", "view", "-c", "-F", "2308", bam_path, ref_name]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            mapped_reads[ref_name] = int(result.stdout.strip())
        except Exception as e:
            print(f"Warning: samtools view failed for {ref_name} in {bam_path}: {e}")
            mapped_reads[ref_name] = 0

    # Run samtools depth
    try:
        cmd = ["samtools", "depth", "-aa", "-J", bam_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                ref_name = parts[0]
                pos = int(parts[1]) - 1 # 0-indexed
                depth = int(parts[2])
                if ref_name in depths:
                    ref_len = ref_lengths[ref_name]
                    if 0 <= pos < ref_len:
                        depths[ref_name][pos] = depth
    except Exception as e:
        print(f"Warning: samtools depth failed for {bam_path}: {e}")
        pass

    # Compute statistics and profile for each reference
    for ref_name, ref_len in ref_lengths.items():
        ref_depths = depths[ref_name]
        mean_depth = sum(ref_depths) / ref_len if ref_len > 0 else 0
        cov_10x = sum(1 for d in ref_depths if d >= 10) / ref_len if ref_len > 0 else 0
        cov_30x = sum(1 for d in ref_depths if d >= 30) / ref_len if ref_len > 0 else 0
        cov_min_cov = sum(1 for d in ref_depths if d >= min_cov) / ref_len if ref_len > 0 else 0
        
        bin_size = max(1, ref_len // num_bins)
        downsampled_depths = []
        for i in range(0, ref_len, bin_size):
            chunk = ref_depths[i : i + bin_size]
            avg_depth = round(sum(chunk) / len(chunk)) if chunk else 0
            downsampled_depths.append({
                "pos": i + 1,
                "depth": avg_depth
            })
            
        results[ref_name] = {
            "mean_depth": round(mean_depth, 2),
            "coverage_breadth_10x": round(cov_10x * 100, 2),
            "coverage_breadth_30x": round(cov_30x * 100, 2),
            "coverage_breadth_min_cov": round(cov_min_cov * 100, 2),
            "mapped_reads": mapped_reads[ref_name],
            "coverage_profile": downsampled_depths
        }
    return results, depths

GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

COMPLEMENT = {
    'A':'T', 'T':'A', 'C':'G', 'G':'C',
    'a':'t', 't':'a', 'c':'g', 'g':'c',
    'N':'N', 'n':'n'
}

def find_matching_model(ref_id, gene, models):
    if not models:
        return None
    for model in models:
        # Check reference
        ref_match = model.get("reference")
        if ref_match and ref_match != ref_id:
            continue
            
        # Check strand
        model_strand = model.get("strand")
        if model_strand and model_strand != gene.get("strand"):
            continue
            
        # Check feature name/label/product/ID
        feature_names = model.get("feature_names", [])
        products = model.get("products", [])
        
        name_match = (
            gene.get("name") in feature_names or
            gene.get("id") in feature_names or
            gene.get("gene") in feature_names or
            gene.get("feature_name") in feature_names or
            gene.get("product") in products
        )
        if not name_match:
            continue
            
        # Check translation exception status
        model_status = model.get("translation_exception_status")
        if model_status and model_status != gene.get("translation_exception_status"):
            continue
            
        return model
    return None

def run_curated_aa_annotation(pos, ref, alt, model, ref_seq):
    if not ref_seq:
        return {"aa_change": "AA annotation unavailable", "status": "curated_model_reference_mismatch"}
        
    strand = model.get("strand", "+")
    
    # 1. Reconstruct the curated coding sequence and trace genomic positions
    cds_bases = []
    cds_positions = []
    prev_end = None
    
    for intv in model.get("intervals", []):
        start = intv.get("start")
        end = intv.get("end")
        phase = intv.get("phase", 0)
        frame_offset = intv.get("frame_offset", 0)
        
        # Check boundary/coordinate validity against reference length
        if start < 1 or end > len(ref_seq) or start > end:
            return {"aa_change": "AA annotation unavailable", "status": "curated_model_reference_mismatch"}
            
        # Validate phase and frame_offset
        if frame_offset not in [-1, 0, 1]:
            return {"aa_change": "AA annotation unavailable", "status": "curated_model_boundary_unsupported"}
            
        if strand == "+":
            intv_pos = list(range(start, end + 1))
            intv_seq = list(ref_seq[start - 1 : end])
            # If interval does not repeat the overlap coordinate but specifies frame_offset == -1
            # (e.g. adjacent intervals start == prev_end + 1 with programmed -1 frameshift):
            if prev_end is not None and frame_offset == -1 and start == prev_end + 1:
                intv_pos = [prev_end] + intv_pos
                intv_seq = [ref_seq[prev_end - 1]] + intv_seq
        else:
            intv_pos = list(range(end, start - 1, -1))
            raw_seq = ref_seq[start - 1 : end]
            intv_seq = [COMPLEMENT.get(b, b) for b in reversed(raw_seq)]
            if prev_end is not None and frame_offset == -1 and start == prev_end + 1:
                intv_pos = [prev_end] + intv_pos
                intv_seq = [COMPLEMENT.get(ref_seq[prev_end - 1], "N")] + intv_seq
                
        prev_end = end
            
        # Apply phase/trim (number of bases to skip from start of interval)
        if phase < 0 or phase > len(intv_seq):
            return {"aa_change": "AA annotation unavailable", "status": "curated_model_boundary_unsupported"}
            
        intv_seq = intv_seq[phase:]
        intv_pos = intv_pos[phase:]
        
        cds_bases.extend(intv_seq)
        cds_positions.extend(intv_pos)
        
    cds_sequence = "".join(cds_bases)
    
    # 2. Map the genomic variant position to curated CDS coordinates
    matching_indices = [i for i, p in enumerate(cds_positions) if p == pos]
    if not matching_indices:
        return {"aa_change": "AA annotation unavailable", "status": "curated_model_position_out_of_bounds"}
        
    # 3. Perform translation for each matching index
    aa_changes = []
    mut_base = alt if strand == "+" else COMPLEMENT.get(alt, alt)
    
    for cds_idx in matching_indices:
        codon_start_idx = (cds_idx // 3) * 3
        codon_offset = cds_idx % 3
        
        # Codon incomplete or out of bounds checks
        if codon_start_idx < 0 or codon_start_idx + 3 > len(cds_sequence):
            return {"aa_change": "AA annotation unavailable", "status": "curated_model_codon_incomplete"}
            
        wt_codon = cds_sequence[codon_start_idx : codon_start_idx + 3]
        if len(wt_codon) != 3:
            return {"aa_change": "AA annotation unavailable", "status": "curated_model_codon_incomplete"}
            
        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = mut_base
        mut_codon = "".join(mut_codon_list)
        
        wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
        mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
        aa_pos = (cds_idx // 3) + 1
        
        if wt_aa == "?" or mut_aa == "?":
            return {"aa_change": "AA annotation unavailable", "status": "curated_model_boundary_unsupported"}
            
        if wt_aa == mut_aa:
            aa_changes.append(f"{wt_aa}{aa_pos}{mut_aa} (syn)")
        else:
            aa_changes.append(f"{wt_aa}{aa_pos}{mut_aa}")
            
    if not aa_changes:
        return {"aa_change": "AA annotation unavailable", "status": "curated_model_boundary_unsupported"}
        
    return {"aa_change": "/".join(aa_changes), "status": None}

def get_aa_change(pos, ref, alt, gene, ref_seq):
    """Calculate codon-level amino acid change for SNPs in coding genes."""
    if not ref_seq or len(ref) != 1 or len(alt) != 1:
        return {"aa_change": "Indel" if len(ref) != len(alt) else "", "status": None}
        
    if gene.get("translation_exception_status"):
        # Lazy load curated models
        models = load_curated_models_cached()
        matching_model = find_matching_model(gene.get("ref_id"), gene, models)
        if not matching_model:
            # If no curated model matches, fail closed with original status
            return {"aa_change": "AA annotation unavailable", "status": gene["translation_exception_status"]}
            
        # If curated model matches, run curated annotation
        return run_curated_aa_annotation(pos, ref, alt, matching_model, ref_seq)
        
    if gene.get("phase", 0) == -1 or gene.get("multi_part", False):
        return {"aa_change": "AA annotation unavailable", "status": None}
        
    g_start = gene["start"]
    g_end = gene["end"]
    strand = gene.get("strand", "+")
    phase = gene.get("phase", 0)
    
    if strand == "+":
        if not (g_start <= pos <= g_end):
            return {"aa_change": "", "status": None}
        rel_offset = (pos - g_start) + phase
        codon_index = rel_offset // 3
        codon_offset = rel_offset % 3
        codon_start = (g_start - 1 - phase) + codon_index * 3
        
        if codon_start < 0 or codon_start + 3 > len(ref_seq):
            return {"aa_change": "", "status": None}
            
        wt_codon = ref_seq[codon_start : codon_start + 3]
        if len(wt_codon) != 3:
            return {"aa_change": "", "status": None}
            
        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = alt
        mut_codon = "".join(mut_codon_list)
        
        wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
        mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
        aa_pos = codon_index + 1
        
        if wt_aa == "?" or mut_aa == "?":
            return {"aa_change": "", "status": None}
        if wt_aa == mut_aa:
            return {"aa_change": f"{wt_aa}{aa_pos}{mut_aa} (syn)", "status": None}
        else:
            return {"aa_change": f"{wt_aa}{aa_pos}{mut_aa}", "status": None}
            
    else: # strand == "-"
        if not (g_start <= pos <= g_end):
            return {"aa_change": "", "status": None}
        rel_offset = (g_end - pos) + phase
        codon_index = rel_offset // 3
        codon_offset = rel_offset % 3
        
        idx1 = (g_end - 1 + phase) - codon_index * 3
        idx3 = idx1 - 2
        
        if idx3 < 0 or idx1 + 1 > len(ref_seq):
            return {"aa_change": "", "status": None}
            
        genomic_codon = ref_seq[idx3 : idx1 + 1]
        if len(genomic_codon) != 3:
            return {"aa_change": "", "status": None}
            
        wt_codon = "".join(COMPLEMENT.get(b, b) for b in reversed(genomic_codon))
        
        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = COMPLEMENT.get(alt, alt)
        mut_codon = "".join(mut_codon_list)
        
        wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
        mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
        aa_pos = codon_index + 1
        
        if wt_aa == "?" or mut_aa == "?":
            return {"aa_change": "", "status": None}
        if wt_aa == mut_aa:
            return {"aa_change": f"{wt_aa}{aa_pos}{mut_aa} (syn)", "status": None}
        else:
            return {"aa_change": f"{wt_aa}{aa_pos}{mut_aa}", "status": None}

def parse_vcf(vcf_path, genes_annotation, ref_seqs=None, min_freq=0.2):
    """Parse VCF file for variants and intersect with GFF genes to find effects, grouping by reference."""
    variants = []
    if not os.path.exists(vcf_path):
        return variants
        
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
                
            ref_chrom = parts[0]
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]
            qual = float(parts[5]) if parts[5] != "." else 30.0
            info_str = parts[7]
            
            # Parse INFO fields
            info = {}
            for item in info_str.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info[k] = v
                else:
                    info[item] = True
            
            # Estimate Allele Frequency (AF) and Depth (DP)
            dp = int(info.get("DP", 100))
            af = float(info.get("AF", 1.0))
            if "DP4" in info:
                # Common in samtools/bcftools
                dp4 = [int(x) for x in info["DP4"].split(",")]
                ref_cnt = dp4[0] + dp4[1]
                alt_cnt = dp4[2] + dp4[3]
                dp = ref_cnt + alt_cnt
                af = alt_cnt / dp if dp > 0 else 1.0
                
            # Strict allele frequency filter
            if af < min_freq:
                continue
                
            # Intersect with GFF genes (supporting overlapping reading frames)
            ref_genes = genes_annotation.get(ref_chrom, [])
            matching_genes = [g for g in ref_genes if g["type"] in ["gene", "CDS"] and g["start"] <= pos <= g["end"]]
            
            if not matching_genes:
                effect_gene = "Intergenic"
                effect_type = "Non-coding"
                aa_change = ""
                status = None
            elif len(matching_genes) == 1:
                gene = matching_genes[0]
                effect_gene = gene["name"]
                effect_type = "Coding"
                if ref_seqs:
                    res = get_aa_change(pos, ref, alt, gene, ref_seqs.get(ref_chrom, ""))
                    aa_change = res["aa_change"]
                    status = res["status"]
                else:
                    aa_change = ""
                    status = None
            else:
                effect_type = "Coding"
                gene_names = list(dict.fromkeys(g["name"] for g in matching_genes))
                effect_gene = "/".join(gene_names)
                aa_changes = []
                statuses = []
                if ref_seqs:
                    for gene in matching_genes:
                        res = get_aa_change(pos, ref, alt, gene, ref_seqs.get(ref_chrom, ""))
                        if res["aa_change"] and res["aa_change"] not in aa_changes:
                            aa_changes.append(res["aa_change"])
                        if res["status"] and res["status"] not in statuses:
                            statuses.append(res["status"])
                aa_change = "/".join(aa_changes) if aa_changes else ""
                status = statuses[0] if statuses else None
                    
            var_dict = {
                "reference": ref_chrom,
                "position": pos,
                "ref": ref,
                "alt": alt,
                "quality": round(qual, 1),
                "depth": dp,
                "frequency": round(af, 3),
                "gene": effect_gene,
                "type": effect_type,
                "aa_change": aa_change
            }
            if status:
                var_dict["status"] = status
            variants.append(var_dict)
            
    return variants

def parse_consensus(fasta_path, sample_id=None):
    """Analyze consensus FASTA sequence for all references."""
    results = {}
    if not os.path.exists(fasta_path):
        return results
        
    consensuses = {}
    current_name = None
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if current_name:
                    consensuses[current_name] = "".join(current_seq)
                current_name = line.strip().split()[0][1:]
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_name:
            consensuses[current_name] = "".join(current_seq)
            
    # Calculate stats for each consensus contig
    for name, seq in consensuses.items():
        ref_name = name
        if sample_id and name.startswith(sample_id + "_"):
            ref_name = name[len(sample_id) + 1:]
            
        length = len(seq)
        g_count = seq.count('G')
        c_count = seq.count('C')
        n_count = seq.count('N')
        
        gc_rate = (g_count + c_count) / (length - n_count) if (length - n_count) > 0 else 0.0
        n_rate = n_count / length if length > 0 else 0.0
        
        results[ref_name] = {
            "length": length,
            "gc_content": round(gc_rate * 100, 2),
            "n_content": round(n_rate * 100, 2),
            "sequence": seq
        }
    return results

def parse_metadata(metadata_path, processed_sample_ids):
    """Parse sample metadata CSV and validate sample IDs."""
    metadata_map = {}
    if not metadata_path or not os.path.exists(metadata_path):
        return metadata_map

    import csv
    with open(metadata_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "sample_id" not in reader.fieldnames:
            print("Error: Metadata CSV must contain a 'sample_id' column.", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            s_id = row.get("sample_id")
            if not s_id:
                continue
            meta_entry = {k: v for k, v in row.items() if k != "sample_id"}
            metadata_map[s_id] = meta_entry

    processed_set = set(processed_sample_ids)
    metadata_set = set(metadata_map.keys())

    unmatched_metadata = metadata_set - processed_set
    if unmatched_metadata:
        print(f"Warning: The following sample IDs in the metadata CSV were not processed in the pipeline run: {', '.join(unmatched_metadata)}", file=sys.stderr)

    unmatched_processed = processed_set - metadata_set
    if unmatched_processed:
        print(f"Warning: The following processed sample IDs had no matching metadata entries in the CSV: {', '.join(unmatched_processed)}", file=sys.stderr)

    return metadata_map

def parse_and_calculate_amplicons(primer_bed_path, amplicons_bed_path, depths, min_cov):
    """Parse amplicon boundaries (directly or via primer BED pairing) and calculate coverage statistics."""
    amplicons = []
    
    if amplicons_bed_path and os.path.exists(amplicons_bed_path):
        with open(amplicons_bed_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    amplicons.append({
                        "chrom": parts[0],
                        "start": int(parts[1]),
                        "end": int(parts[2]),
                        "name": parts[3]
                    })
        return calculate_amplicon_coverage(amplicons, depths, min_cov)
        
    if primer_bed_path and os.path.exists(primer_bed_path):
        import re
        suffix_pat = re.compile(r'[-_](LEFT|RIGHT|FORWARD|REVERSE|F|R|L|R|left|right|forward|reverse|f|r)(_alt\d*)?$', re.IGNORECASE)
        primers = defaultdict(list)
        
        with open(primer_bed_path, 'r') as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3]
                    prefix = suffix_pat.sub('', name)
                    primers[(chrom, prefix)].append((start, end))
                    
        for (chrom, prefix), p_list in primers.items():
            if len(p_list) < 2:
                continue
            start = min(p[0] for p in p_list)
            end = max(p[1] for p in p_list)
            amplicons.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "name": prefix
            })
            
        if not amplicons:
            print("Error: Could not reliably pair primers into amplicons from the primer BED file. Please check naming conventions (e.g. using '_LEFT' / '_RIGHT' or '_F' / '_R' suffixes) or supply an explicit amplicon coordinates file using --amplicons-bed.", file=sys.stderr)
            sys.exit(1)
            
        return calculate_amplicon_coverage(amplicons, depths, min_cov)
        
    return []

def calculate_amplicon_coverage(amplicons, depths, min_cov):
    """Helper to calculate per-amplicon coverage stats using depths dictionary."""
    amplicon_metrics = []
    for amp in amplicons:
        chrom = amp["chrom"]
        start = amp["start"]
        end = amp["end"]
        name = amp["name"]
        
        chrom_depths = depths.get(chrom, [])
        if not chrom_depths:
            continue
            
        amp_depths = chrom_depths[start : min(end, len(chrom_depths))]
        if not amp_depths:
            continue
            
        mean_depth = sum(amp_depths) / len(amp_depths)
        min_depth = min(amp_depths)
        
        sorted_depths = sorted(amp_depths)
        median_depth = sorted_depths[len(sorted_depths) // 2]
        
        bases_above = sum(1 for d in amp_depths if d >= min_cov)
        pct_above = (bases_above / len(amp_depths)) * 100
        
        status = "pass" if pct_above >= 90.0 or mean_depth >= min_cov else "fail"
        
        amplicon_metrics.append({
            "amplicon_id": name,
            "reference": chrom,
            "start": start + 1,
            "end": end,
            "mean_depth": round(mean_depth, 2),
            "min_depth": min_depth,
            "median_depth": median_depth,
            "percent_bases_above_threshold": round(pct_above, 2),
            "status": status
        })
    return amplicon_metrics

def calculate_distance_matrices(samples_data, ref_info, min_cov):
    """Calculate coordinate-aligned SNP-only pairwise distance matrices for each reference strain over the entire genome."""
    distance_matrices = {}
    samples_list = list(samples_data.values())
    num_samples = len(samples_list)
    distance_matrix_order = [s["sample_id"] for s in samples_list]
    distance_matrices_by_name = {}
    
    for ref in ref_info:
        ref_name = ref["name"]
        ref_len = ref["length"]
        ref_seq = ref["sequence"]
        
        # Initialize N x N matrix with 0.0
        matrix = [[0.0] * num_samples for _ in range(num_samples)]
        
        # 1. Compare every pair of samples base-by-base across the ENTIRE reference genome
        for i in range(num_samples):
            sA = samples_list[i]
            depths_A = sA.get("full_depths", {}).get(ref_name, [0] * ref_len)
            # Map of positions (1-indexed) to their alternative alleles (major consensus SNPs only)
            vars_A = {
                v["position"]: v["alt"] 
                for v in sA.get("variants", []) 
                if v["reference"] == ref_name and len(v["ref"]) == 1 and len(v["alt"]) == 1 and v.get("frequency", 1.0) >= 0.5
            }
            
            # Track positions deleted in sample A to exclude from SNP comparison
            del_positions_A = set()
            for v in sA.get("variants", []):
                if v["reference"] == ref_name and len(v["ref"]) > len(v["alt"]) and v.get("frequency", 1.0) >= 0.5:
                    start_del = v["position"] + len(v["alt"])
                    end_del = v["position"] + len(v["ref"]) - 1
                    del_positions_A.update(range(start_del, end_del + 1))
            
            for j in range(i + 1, num_samples):
                sB = samples_list[j]
                depths_B = sB.get("full_depths", {}).get(ref_name, [0] * ref_len)
                vars_B = {
                    v["position"]: v["alt"] 
                    for v in sB.get("variants", []) 
                    if v["reference"] == ref_name and len(v["ref"]) == 1 and len(v["alt"]) == 1 and v.get("frequency", 1.0) >= 0.5
                }
                
                del_positions_B = set()
                for v in sB.get("variants", []):
                    if v["reference"] == ref_name and len(v["ref"]) > len(v["alt"]) and v.get("frequency", 1.0) >= 0.5:
                        start_del = v["position"] + len(v["alt"])
                        end_del = v["position"] + len(v["ref"]) - 1
                        del_positions_B.update(range(start_del, end_del + 1))
                
                diffs = 0
                compared = 0
                
                for pos_idx in range(ref_len):
                    pos = pos_idx + 1
                    
                    # Skip positions deleted in either sample so deletions are not falsely counted as reference matches
                    if pos in del_positions_A or pos in del_positions_B:
                        continue
                        
                    depth_A = depths_A[pos_idx]
                    depth_B = depths_B[pos_idx]
                    
                    # Only compare if both samples have sufficient coverage depth
                    if depth_A >= min_cov and depth_B >= min_cov:
                        compared += 1
                        base_A = vars_A.get(pos, ref_seq[pos_idx])
                        base_B = vars_B.get(pos, ref_seq[pos_idx])
                        if base_A.upper() != base_B.upper():
                            diffs += 1
                            
                dist = diffs / compared if compared > 0 else 1.0
                matrix[i][j] = dist
                matrix[j][i] = dist
                
        distance_matrices[ref_name] = matrix
        
        # Construct name-keyed distance map for ref_name
        ref_by_name = {}
        for i, sA_id in enumerate(distance_matrix_order):
            ref_by_name[sA_id] = {}
            for j, sB_id in enumerate(distance_matrix_order):
                ref_by_name[sA_id][sB_id] = matrix[i][j]
        distance_matrices_by_name[ref_name] = ref_by_name
        
    return distance_matrices, distance_matrix_order, distance_matrices_by_name

def validate_report_schema(report):
    """
    Validates that the output report dictionary matches the required JSON Schema
    structure and data types for the Viromics Dashboard, ensuring zero-dependency safety.
    """
    # 1. Top-Level Structure
    if not isinstance(report, dict):
        raise ValueError("Schema Error: Report must be a JSON object dictionary.")
        
    required_top = ["run_id", "min_coverage_threshold", "references", "genes", "distance_matrices", "distance_matrix_order", "distance_matrices_by_name", "samples"]
    for key in required_top:
        if key not in report:
            raise ValueError(f"Schema Error: Report missing required top-level key '{key}'.")
            
    if not isinstance(report["run_id"], str):
        raise ValueError("Schema Error: 'run_id' must be a string.")
    if not isinstance(report["min_coverage_threshold"], int):
        raise ValueError("Schema Error: 'min_coverage_threshold' must be an integer.")
    if not isinstance(report["references"], list):
        raise ValueError("Schema Error: 'references' must be a list.")
    if not isinstance(report["genes"], dict):
        raise ValueError("Schema Error: 'genes' must be an object dictionary.")
    if not isinstance(report["distance_matrices"], dict):
        raise ValueError("Schema Error: 'distance_matrices' must be an object dictionary.")
    if not isinstance(report["distance_matrix_order"], list):
        raise ValueError("Schema Error: 'distance_matrix_order' must be a list.")
    for idx, val in enumerate(report["distance_matrix_order"]):
        if not isinstance(val, str):
            raise ValueError(f"Schema Error: 'distance_matrix_order' entry at index {idx} must be a string.")
    if not isinstance(report["distance_matrices_by_name"], dict):
        raise ValueError("Schema Error: 'distance_matrices_by_name' must be an object dictionary.")
    for ref_name, matrix_dict in report["distance_matrices_by_name"].items():
        if not isinstance(matrix_dict, dict):
            raise ValueError(f"Schema Error: 'distance_matrices_by_name' for ref '{ref_name}' must be a dictionary.")
        for sA_id, sB_dict in matrix_dict.items():
            if not isinstance(sB_dict, dict):
                raise ValueError(f"Schema Error: 'distance_matrices_by_name' entry for ref '{ref_name}', sample '{sA_id}' must be a dictionary.")
            for sB_id, dist in sB_dict.items():
                if not isinstance(dist, (int, float)):
                    raise ValueError(f"Schema Error: 'distance_matrices_by_name' distance for ref '{ref_name}', samples '{sA_id}' vs '{sB_id}' must be a number.")
    if not isinstance(report["samples"], list):
        raise ValueError("Schema Error: 'samples' must be a list.")
        
    # 2. References list
    for idx, ref in enumerate(report["references"]):
        if not isinstance(ref, dict):
            raise ValueError(f"Schema Error: Reference entry at index {idx} is not a dictionary.")
        for key in ["name", "length", "sequence"]:
            if key not in ref:
                raise ValueError(f"Schema Error: Reference at index {idx} missing key '{key}'.")
        if not isinstance(ref["name"], str):
            raise ValueError(f"Schema Error: Reference at index {idx} 'name' must be a string.")
        if not isinstance(ref["length"], int):
            raise ValueError(f"Schema Error: Reference at index {idx} 'length' must be an integer.")
        if not isinstance(ref["sequence"], str):
            raise ValueError(f"Schema Error: Reference at index {idx} 'sequence' must be a string.")
            
    # 3. Samples list
    for idx, s in enumerate(report["samples"]):
        if not isinstance(s, dict):
            raise ValueError(f"Schema Error: Sample at index {idx} is not a dictionary.")
        required_sample = ["sample_id", "qc", "coverage", "variants", "consensus"]
        for key in required_sample:
            if key not in s:
                raise ValueError(f"Schema Error: Sample at index {idx} (ID: {s.get('sample_id', 'unknown')}) missing required key '{key}'.")
                
        s_id = s["sample_id"]
        if not isinstance(s_id, str):
            raise ValueError(f"Schema Error: Sample at index {idx} 'sample_id' must be a string.")
            
        # Sample QC validation
        qc = s["qc"]
        if not isinstance(qc, dict):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'qc' must be a dictionary.")
        required_qc = ["raw_reads", "raw_bases", "trimmed_reads", "trimmed_bases", "q30_rate_raw", "q30_rate_trimmed", "duplication_rate"]
        for key in required_qc:
            if key not in qc:
                raise ValueError(f"Schema Error: Sample '{s_id}' QC missing key '{key}'.")
            if not isinstance(qc[key], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' QC key '{key}' must be a number.")
                
        # Sample Coverage validation
        coverage = s["coverage"]
        if not isinstance(coverage, dict):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'coverage' must be a dictionary.")
        for ref_name, cov_data in coverage.items():
            if not isinstance(cov_data, dict):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage data for ref '{ref_name}' must be a dictionary.")
            required_cov = ["mean_depth", "coverage_breadth_10x", "coverage_breadth_30x", "coverage_breadth_min_cov", "mapped_reads", "coverage_profile"]
            for key in required_cov:
                if key not in cov_data:
                    raise ValueError(f"Schema Error: Sample '{s_id}' coverage for ref '{ref_name}' missing key '{key}'.")
            if not isinstance(cov_data["mean_depth"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage mean_depth for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["coverage_breadth_10x"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage breadth_10x for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["coverage_breadth_30x"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage breadth_30x for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["coverage_breadth_min_cov"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage breadth_min_cov for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["mapped_reads"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage mapped_reads for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["coverage_profile"], list):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage_profile for ref '{ref_name}' must be a list.")
                
        # Sample Consensus validation
        consensus = s["consensus"]
        if not isinstance(consensus, dict):
            raise ValueError(f"Schema Error: Sample '{s_id}' 'consensus' must be a dictionary.")
        for ref_name, cons_data in consensus.items():
            if not isinstance(cons_data, dict):
                raise ValueError(f"Schema Error: Sample '{s_id}' consensus data for ref '{ref_name}' must be a dictionary.")
            required_cons = ["length", "gc_content", "n_content", "sequence"]
            for key in required_cons:
                if key not in cons_data:
                    raise ValueError(f"Schema Error: Sample '{s_id}' consensus for ref '{ref_name}' missing key '{key}'.")
            if not isinstance(cons_data["length"], int):
                raise ValueError(f"Schema Error: Sample '{s_id}' consensus length for ref '{ref_name}' must be an integer.")
            if not isinstance(cons_data["gc_content"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' consensus gc_content for ref '{ref_name}' must be a number.")
            if not isinstance(cons_data["n_content"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' consensus n_content for ref '{ref_name}' must be a number.")
            if not isinstance(cons_data["sequence"], str):
                raise ValueError(f"Schema Error: Sample '{s_id}' consensus sequence for ref '{ref_name}' must be a string.")

def main():
    args = parse_args()
    
    # Parse reference genome (can contain multiple segments/strains)
    ref_info = parse_reference(args.reference)
    ref_lengths = {ref["name"]: ref["length"] for ref in ref_info}
    ref_seqs = {ref["name"]: ref["sequence"] for ref in ref_info}
    
    # Parse GFF annotations (grouped by sequence ID)
    genes = parse_gff(args.gff)
    
    # Split input paths
    fastp_files = args.fastp.split(",")
    bam_files = args.bams.split(",")
    vcf_files = args.vcfs.split(",")
    consensus_files = args.consensuses.split(",")
    
    # Map sample IDs
    samples_data = {}
    
    for fastp_f in fastp_files:
        sample_id = os.path.basename(fastp_f).replace("_fastp.json", "")
        samples_data[sample_id] = {
            "sample_id": sample_id,
            "qc": parse_fastp(fastp_f)
        }
        
    for bam_f in bam_files:
        sample_id = os.path.basename(bam_f).replace("_trimmed.bam", "").replace("_raw.bam", "")
        if sample_id in samples_data:
            cov, full_depths = calculate_coverage(bam_f, ref_lengths, min_cov=args.min_coverage)
            samples_data[sample_id]["coverage"] = cov
            samples_data[sample_id]["full_depths"] = full_depths
            
    for vcf_f in vcf_files:
        sample_id = os.path.basename(vcf_f).replace(".vcf", "")
        if sample_id in samples_data:
            variants = parse_vcf(vcf_f, genes, ref_seqs, args.min_freq)
            samples_data[sample_id]["variants"] = variants
            
    for cons_f in consensus_files:
        sample_id = os.path.basename(cons_f).replace("_consensus.fasta", "")
        if sample_id in samples_data:
            cons_info = parse_consensus(cons_f, sample_id)
            samples_data[sample_id]["consensus"] = cons_info

    # Parse metadata CSV
    metadata_map = parse_metadata(args.metadata, list(samples_data.keys()))
    
    # Calculate amplicon metrics and embed metadata
    for sample_id, s_data in samples_data.items():
        s_data["metadata"] = metadata_map.get(sample_id, {})
        s_data["amplicons"] = []
        if (args.primers or args.amplicons_bed) and "full_depths" in s_data:
            s_data["amplicons"] = parse_and_calculate_amplicons(
                args.primers, args.amplicons_bed, s_data["full_depths"], args.min_coverage
            )

    # Calculate pairwise distance matrices on pipeline side (SNP-only, frame-shift-proof)
    distance_matrices, distance_matrix_order, distance_matrices_by_name = calculate_distance_matrices(samples_data, ref_info, args.min_coverage)
    
    # Remove bulky full_depths maps to keep JSON report file size compact
    for s in samples_data.values():
        if "full_depths" in s:
            del s["full_depths"]
            
    # Aggregate into final report structure
    report = {
        "run_id": args.run_id,
        "min_coverage_threshold": args.min_coverage,
        "references": ref_info,
        "genes": genes,
        "distance_matrices": distance_matrices,
        "distance_matrix_order": distance_matrix_order,
        "distance_matrices_by_name": distance_matrices_by_name,
        "samples": list(samples_data.values())
    }
    
    # Validate report against schema definitions
    validate_report_schema(report)
    
    # Write JSON report
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
        
    print(f"Successfully generated Viromics Run Report: {args.output}")

if __name__ == "__main__":
    main()
