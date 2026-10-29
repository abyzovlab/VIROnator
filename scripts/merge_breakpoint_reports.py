#!/usr/bin/env python3
"""
merge_breakpoint_reports.py
Master Consolidation Script for VIROnator Module 10 (Discordant Reads & Integration Breakpoints)

Aggregates individual per-sample breakpoint TSVs (<SAMPLE>_integration_breakpoints.tsv)
from the vironator output directory into a single master cohort report:
cohort_breakpoints_master.tsv
"""

import os
import sys
import glob
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Consolidate per-sample integration breakpoint TSVs into a master cohort report.")
    parser.add_argument(
        "--vironator-dir",
        default=os.environ.get("VIRONATOR_DIR", "/mnt/disks/staff/vironator"),
        help="Base directory containing per-sample vironator outputs"
    )
    parser.add_argument(
        "--out-file",
        default="cohort_breakpoints_master.tsv",
        help="Output path for master breakpoint report TSV"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header = [
        "sample_id", "phase", "project", "human_chr", "human_coord_window", "human_breakpoint_est",
        "viral_accession", "viral_name_sanitized", "viral_coord_window", "viral_breakpoint_est",
        "orientation", "discordant_read_count", "split_read_count", "total_evidence_reads",
        "confidence", "intragenic_element", "upstream_element", "upstream_distance_bp",
        "downstream_element", "downstream_distance_bp"
    ]

    pattern = os.path.join(args.vironator_dir, "**", "*_integration_breakpoints.tsv")
    tsv_files = sorted(glob.glob(pattern, recursive=True))

    print(f"[INFO] Scanning for breakpoint reports in {args.vironator_dir} ...")
    print(f"[INFO] Found {len(tsv_files)} per-sample breakpoint TSV files.")

    seen = set()
    rows = []
    for fpath in tsv_files:
        if os.path.basename(fpath) == os.path.basename(args.out_file):
            continue
        try:
            with open(fpath, "r", encoding="utf-8", errors="replace") as f:
                lines = [line.strip() for line in f if line.strip()]
                if len(lines) >= 2:
                    for line in lines[1:]:
                        parts = line.split("\t")
                        if len(parts) >= len(header):
                            # Deduplication key: sample_id, human_chr, human_bp, viral_acc, viral_bp
                            key = (parts[0], parts[3], parts[5], parts[6], parts[9])
                            if key not in seen:
                                seen.add(key)
                                rows.append(parts)
        except Exception as e:
            print(f"[WARNING] Could not read {fpath}: {e}")

    # Sort rows by sample_id, then total_evidence_reads descending
    try:
        rows.sort(key=lambda r: (r[0], -int(r[13]) if r[13].isdigit() else 0))
    except Exception:
        rows.sort(key=lambda r: r[0])

    with open(args.out_file, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Consolidated {len(rows)} breakpoint records into {args.out_file}")


if __name__ == "__main__":
    main()
