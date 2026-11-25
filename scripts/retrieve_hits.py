#!/usr/bin/env python3
import os
import sys
import json
import argparse
import csv
import re

def sanitize_value(val, delimiter):
    """Sanitize metadata values to be FASTA-header friendly."""
    if val is None:
        return "unknown"
    val = str(val).strip()
    if not val:
        return "unknown"
    # Replace spaces and non-alphanumeric chars (excluding dots, hyphens, and underscores) with delimiter
    sanitized = re.sub(r'[^a-zA-Z0-9._-]', delimiter, val)
    # Collapse multiple delimiters
    sanitized = re.sub(re.escape(delimiter) + r'+', delimiter, sanitized)
    return sanitized.strip(delimiter)

def parse_custom_metadata(csv_path):
    """Parse custom metadata sheet if provided."""
    metadata_map = {}
    if not csv_path or not os.path.exists(csv_path):
        return None
    
    with open(csv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or 'sample_id' not in reader.fieldnames:
            print(f"Error: Custom metadata CSV {csv_path} must contain a 'sample_id' column.", file=sys.stderr)
            return None
        for row in reader:
            sample_id = row.get('sample_id')
            if sample_id:
                metadata_map[sample_id] = row
    return metadata_map

def parse_fasta_file(fasta_path):
    """Simple parser to read FASTA records from file."""
    records = {}
    if not os.path.exists(fasta_path):
        return records
    
    current_header = None
    current_seq = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    records[current_header] = "".join(current_seq)
                current_header = line[1:]
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_header:
            records[current_header] = "".join(current_seq)
    return records

def main():
    parser = argparse.ArgumentParser(
        description="Extract consensus sequences for mapped reference contigs exceeding a coverage breadth threshold."
    )
    parser.add_argument(
        "-r", "--report", required=True,
        help="Path to viromics_run_report.json file."
    )
    parser.add_argument(
        "-o", "--output", default="combined_hits.fasta",
        help="Path to write combined FASTA output (default: combined_hits.fasta)."
    )
    parser.add_argument(
        "-b", "--min-breadth", type=float, default=10.0,
        help="Minimum coverage breadth threshold (0-100%%) to qualify as a hit (default: 10.0%%)."
    )
    parser.add_argument(
        "-m", "--metadata", default=None,
        help="Optional path to custom metadata.csv (overrides report embedded metadata)."
    )
    parser.add_argument(
        "-f", "--fields", default=None,
        help="Comma-separated list of metadata columns to append to headers (e.g. unique_id,host)."
    )
    parser.add_argument(
        "-d", "--delimiter", default="_",
        help="Delimiter to separate appended header metadata fields (default: _)."
    )
    parser.add_argument(
        "--results-dir", default=None,
        help="Root results folder containing sample results directories (defaults to directory of JSON report)."
    )
    parser.add_argument(
        "-j", "--use-json", action="store_true",
        help="Extract sequences directly from JSON report instead of consensus FASTA files."
    )

    args = parser.parse_args()

    # Verify input report
    if not os.path.exists(args.report):
        print(f"Error: Report JSON file not found at {args.report}", file=sys.stderr)
        sys.exit(1)

    # Load report
    try:
        with open(args.report, 'r') as f:
            report = json.load(f)
    except Exception as e:
        print(f"Error: Failed to parse report JSON: {e}", file=sys.stderr)
        sys.exit(1)

    # Verify basic JSON structure
    if "samples" not in report:
        print("Error: Invalid report JSON schema (missing 'samples' array).", file=sys.stderr)
        sys.exit(1)

    # Load custom metadata if specified
    custom_metadata = parse_custom_metadata(args.metadata) if args.metadata else None

    # Deduce results directory
    results_dir = args.results_dir
    if not results_dir:
        results_dir = os.path.dirname(os.path.abspath(args.report))

    # Parse metadata columns to append
    metadata_fields = [f.strip() for f in args.fields.split(',')] if args.fields else []

    extracted_records = []
    skipped_count = 0
    hit_count = 0

    print(f"Processing report '{args.report}'...")
    print(f"Filtering hits with coverage breadth >= {args.min_breadth}%...")

    for sample in report["samples"]:
        sample_id = sample.get("sample_id")
        if not sample_id:
            continue

        coverage = sample.get("coverage", {})
        consensus = sample.get("consensus", {})

        # Determine metadata sheet to retrieve values from
        meta_dict = {}
        if custom_metadata:
            meta_dict = custom_metadata.get(sample_id, {})
        else:
            meta_dict = sample.get("metadata", {})

        # Process each contig mapping
        for ref_name, cov_stats in coverage.items():
            # Get coverage breadth at the user-defined minimum coverage threshold (falling back to 10x if absent)
            breadth = cov_stats.get("coverage_breadth_min_cov", cov_stats.get("coverage_breadth_10x", 0.0))

            if breadth >= args.min_breadth:
                hit_count += 1
                sequence = None

                # Method 1: Retrieve from consensus files
                if not args.use_json:
                    cons_fasta = os.path.join(results_dir, sample_id, "consensus", f"{sample_id}_consensus.fasta")
                    if os.path.exists(cons_fasta):
                        records = parse_fasta_file(cons_fasta)
                        # Find matching record in FASTA (the header matches `{sample_id}_{ref_name}`)
                        # Since headers are stored without '>' in parse_fasta_file keys, match directly
                        target_key = f"{sample_id}_{ref_name}"
                        sequence = records.get(target_key)
                        
                        # Fallback match if naming varies slightly
                        if not sequence:
                            for key, seq in records.items():
                                if ref_name in key:
                                    sequence = seq
                                    break
                    else:
                        print(f"Warning: Consensus FASTA file not found for sample {sample_id} at {cons_fasta}. Falling back to JSON report sequence.")

                # Method 2: Fallback to JSON report sequence if consensus file missing or explicitly requested
                if not sequence:
                    ref_cons = consensus.get(ref_name, {})
                    sequence = ref_cons.get("sequence")

                if not sequence:
                    print(f"Warning: Could not retrieve consensus sequence for sample {sample_id} contig {ref_name}.", file=sys.stderr)
                    continue

                # Build new header
                base_header = f"{sample_id}_{ref_name}"
                header_parts = [base_header]

                # Append metadata fields
                for field in metadata_fields:
                    val = meta_dict.get(field)
                    sanitized_val = sanitize_value(val, args.delimiter)
                    header_parts.append(sanitized_val)

                new_header = args.delimiter.join(header_parts)
                extracted_records.append((new_header, sequence))
            else:
                skipped_count += 1

    # Write output FASTA
    if extracted_records:
        try:
            with open(args.output, 'w') as out_f:
                for header, seq in extracted_records:
                    out_f.write(f">{header}\n")
                    # Wrap sequence lines at 60 characters
                    for idx in range(0, len(seq), 60):
                        out_f.write(f"{seq[idx:idx+60]}\n")
            print(f"Success: Extracted {len(extracted_records)} hit sequences (skipped {skipped_count} non-hits).")
            print(f"Combined hits FASTA written to: {args.output}")
        except Exception as e:
            print(f"Error: Failed to write output FASTA file: {e}", file=sys.stderr)
            sys.exit(1)
    else:
        print("Notice: No sequences matched the minimum coverage breadth threshold. No output file written.")

if __name__ == "__main__":
    main()
