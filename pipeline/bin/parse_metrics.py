#!/usr/bin/env python3
import os
import sys
import json
import argparse
import subprocess
from collections import defaultdict

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
    genes = defaultdict(list)
    if not gff_path or not os.path.exists(gff_path):
        return dict(genes)
        
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
                
                # Parse attributes
                attrs = {}
                for attr in parts[8].split(";"):
                    if "=" in attr:
                        k, v = attr.split("=", 1)
                        attrs[k.strip()] = v.strip()
                
                name = attrs.get("Name", attrs.get("gene", attrs.get("ID", "Unknown")))
                product = attrs.get("product", "").replace("%20", " ")
                
                genes[ref_id].append({
                    "name": name,
                    "type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "product": product
                })
    
    # Sort genes for each reference by start position
    sorted_genes = {}
    for ref_id, g_list in genes.items():
        sorted_genes[ref_id] = sorted(g_list, key=lambda x: x["start"])
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
    
    # Run samtools idxstats to retrieve mapped reads per contig
    try:
        cmd = ["samtools", "idxstats", bam_path]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        for line in result.stdout.strip().split("\n"):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                ref_name = parts[0]
                mapped_cnt = int(parts[2])
                if ref_name in mapped_reads:
                    mapped_reads[ref_name] = mapped_cnt
    except Exception as e:
        print(f"Warning: samtools idxstats failed for {bam_path}: {e}")
        pass

    # Run samtools depth
    try:
        cmd = ["samtools", "depth", "-a", bam_path]
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
        cov_10x = sum(1 for d in ref_depths if d >= min_cov) / ref_len if ref_len > 0 else 0
        cov_30x = sum(1 for d in ref_depths if d >= 30) / ref_len if ref_len > 0 else 0
        
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

def get_aa_change(pos, ref, alt, gene, ref_seq):
    """Calculate codon-level amino acid change for SNPs in coding genes."""
    if not ref_seq or len(ref) != 1 or len(alt) != 1:
        return "Indel" if len(ref) != len(alt) else ""
        
    g_start = gene["start"]
    g_end = gene["end"]
    strand = gene.get("strand", "+")
    
    if strand == "+":
        if not (g_start <= pos <= g_end):
            return ""
        rel_offset = pos - g_start
        codon_index = rel_offset // 3
        codon_offset = rel_offset % 3
        codon_start = (g_start - 1) + codon_index * 3
        
        if codon_start < 0 or codon_start + 3 > len(ref_seq):
            return ""
            
        wt_codon = ref_seq[codon_start : codon_start + 3]
        if len(wt_codon) != 3:
            return ""
            
        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = alt
        mut_codon = "".join(mut_codon_list)
        
        wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
        mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
        aa_pos = codon_index + 1
        
        if wt_aa == "?" or mut_aa == "?":
            return ""
        if wt_aa == mut_aa:
            return f"{wt_aa}{aa_pos}{mut_aa} (syn)"
        else:
            return f"{wt_aa}{aa_pos}{mut_aa}"
            
    else: # strand == "-"
        if not (g_start <= pos <= g_end):
            return ""
        rel_offset = g_end - pos
        codon_index = rel_offset // 3
        codon_offset = rel_offset % 3
        
        idx1 = g_end - 1 - codon_index * 3
        idx3 = idx1 - 2
        
        if idx3 < 0 or idx1 + 1 > len(ref_seq):
            return ""
            
        genomic_codon = ref_seq[idx3 : idx1 + 1]
        if len(genomic_codon) != 3:
            return ""
            
        wt_codon = "".join(COMPLEMENT.get(b, b) for b in reversed(genomic_codon))
        
        mut_codon_list = list(wt_codon)
        mut_codon_list[codon_offset] = COMPLEMENT.get(alt, alt)
        mut_codon = "".join(mut_codon_list)
        
        wt_aa = GENETIC_CODE.get(wt_codon.upper(), "?")
        mut_aa = GENETIC_CODE.get(mut_codon.upper(), "?")
        aa_pos = codon_index + 1
        
        if wt_aa == "?" or mut_aa == "?":
            return ""
        if wt_aa == mut_aa:
            return f"{wt_aa}{aa_pos}{mut_aa} (syn)"
        else:
            return f"{wt_aa}{aa_pos}{mut_aa}"

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
                
            # Intersect with GFF genes
            effect_gene = "Intergenic"
            effect_type = "Non-coding"
            aa_change = ""
            ref_genes = genes_annotation.get(ref_chrom, [])
            for gene in ref_genes:
                if gene["type"] in ["gene", "CDS"] and gene["start"] <= pos <= gene["end"]:
                    effect_gene = gene["name"]
                    effect_type = "Coding"
                    if ref_seqs:
                        aa_change = get_aa_change(pos, ref, alt, gene, ref_seqs.get(ref_chrom, ""))
                    break
                    
            variants.append({
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
            })
            
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
        suffix_pat = re.compile(r'_(LEFT|RIGHT|FORWARD|REVERSE|F|R|L|R|left|right|forward|reverse|f|r)(_alt\d*)?$', re.IGNORECASE)
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
            # Map of positions (1-indexed) to their alternative alleles (SNPs only)
            vars_A = {v["position"]: v["alt"] for v in sA.get("variants", []) if v["reference"] == ref_name and len(v["ref"]) == 1 and len(v["alt"]) == 1}
            
            for j in range(i + 1, num_samples):
                sB = samples_list[j]
                depths_B = sB.get("full_depths", {}).get(ref_name, [0] * ref_len)
                vars_B = {v["position"]: v["alt"] for v in sB.get("variants", []) if v["reference"] == ref_name and len(v["ref"]) == 1 and len(v["alt"]) == 1}
                
                diffs = 0
                compared = 0
                
                for pos_idx in range(ref_len):
                    depth_A = depths_A[pos_idx]
                    depth_B = depths_B[pos_idx]
                    
                    # Only compare if both samples have sufficient coverage depth
                    if depth_A >= min_cov and depth_B >= min_cov:
                        compared += 1
                        pos = pos_idx + 1
                        base_A = vars_A.get(pos, ref_seq[pos_idx])
                        base_B = vars_B.get(pos, ref_seq[pos_idx])
                        if base_A.upper() != base_B.upper():
                            diffs += 1
                            
                dist = diffs / compared if compared > 0 else 0.0
                matrix[i][j] = dist
                matrix[j][i] = dist
                
        distance_matrices[ref_name] = matrix
        
    return distance_matrices

def validate_report_schema(report):
    """
    Validates that the output report dictionary matches the required JSON Schema
    structure and data types for the Viromics Dashboard, ensuring zero-dependency safety.
    """
    # 1. Top-Level Structure
    if not isinstance(report, dict):
        raise ValueError("Schema Error: Report must be a JSON object dictionary.")
        
    required_top = ["run_id", "min_coverage_threshold", "references", "genes", "distance_matrices", "samples"]
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
            required_cov = ["mean_depth", "coverage_breadth_10x", "coverage_breadth_30x", "mapped_reads", "coverage_profile"]
            for key in required_cov:
                if key not in cov_data:
                    raise ValueError(f"Schema Error: Sample '{s_id}' coverage for ref '{ref_name}' missing key '{key}'.")
            if not isinstance(cov_data["mean_depth"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage mean_depth for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["coverage_breadth_10x"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage breadth_10x for ref '{ref_name}' must be a number.")
            if not isinstance(cov_data["coverage_breadth_30x"], (int, float)):
                raise ValueError(f"Schema Error: Sample '{s_id}' coverage breadth_30x for ref '{ref_name}' must be a number.")
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
    distance_matrices = calculate_distance_matrices(samples_data, ref_info, args.min_coverage)
    
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
