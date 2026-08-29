#!/usr/bin/env python3
"""
merge_refinement_reports.py
Master Consolidation Script for VIROnator Refinement Module

Aggregates individual per-sample refined report TSVs into a single master cohort report:
cohort_refinement_master.tsv
"""

import os
import sys
import glob
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Consolidate per-sample refined report TSVs into master cohort report.")
    parser.add_argument("--refinement-dir", default="/mnt/disks/staff/SSC_hg38_refinement", help="Base directory containing per-sample refinement outputs")
    parser.add_argument("--out-file", default="cohort_refinement_master.tsv", help="Output path for master refinement report TSV")
    return parser.parse_args()


def main():
    args = parse_args()

    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header = [
        "Sample_ID", "Phase", "Project", "Virus_Accession",
        "Species_TaxID", "Species_Name", "Reference_Unique_Reads",
        "Species_Supportive_Reads", "Total_Refined_Reads",
        "Refined_Copy_Number", "Classification_Status"
    ]

    pattern = os.path.join(args.refinement_dir, "**", "*_refined_report.tsv")
    tsv_files = sorted(glob.glob(pattern, recursive=True))

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
                        if len(parts) >= 11:
                            key = (parts[0], parts[1], parts[2], parts[3], parts[4])
                            if key not in seen:
                                seen.add(key)
                                rows.append(parts)
        except Exception as e:
            print(f"[WARNING] Could not read {fpath}: {e}")

    with open(args.out_file, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Consolidated {len(rows)} refined records from {len(tsv_files)} files into {args.out_file}")


if __name__ == "__main__":
    main()
