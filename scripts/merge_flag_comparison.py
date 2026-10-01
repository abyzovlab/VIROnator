#!/usr/bin/env python3
"""
Automated Consolidation Script for Flag Comparison Module
Merges all individual per-sample flag comparison TSV files into a single master cohort TSV file.
"""

import os
import sys
import glob
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Consolidate Flag Comparison TSV Reports")
    parser.add_argument("--flag-dir", required=True, help="Base directory containing flag comparison sample outputs")
    parser.add_argument("--out-file", required=True, help="Master output TSV report file path")
    return parser.parse_args()


def main():
    args = parse_args()
    
    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header = [
        "phase",
        "project",
        "sample",
        "strategy",
        "original_command_flags_reads",
        "original_command_flags_modified_reads",
        "manual_flags_reads"
    ]

    tsv_pattern = os.path.join(args.flag_dir, "**", "*_flag_comparison.tsv")
    tsv_files = sorted(glob.glob(tsv_pattern, recursive=True))

    seen = set()
    rows = []
    for fpath in tsv_files:
        if os.path.basename(fpath) == os.path.basename(args.out_file):
            continue
        try:
            with open(fpath, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
                if len(lines) >= 2:
                    for line in lines[1:]:
                        parts = line.split("\t")
                        if len(parts) >= 7:
                            # Unique key per (phase, project, sample, strategy)
                            key = (parts[0], parts[1], parts[2], parts[3])
                            if key not in seen:
                                seen.add(key)
                                rows.append(parts)
        except Exception as e:
            print(f"[WARNING] Could not read {fpath}: {e}")

    with open(args.out_file, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Consolidated {len(rows)} unique (phase, project, sample) records from {len(tsv_files)} files into {args.out_file}")


if __name__ == "__main__":
    main()
