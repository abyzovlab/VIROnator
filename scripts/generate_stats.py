#!/usr/bin/env python3
"""
generate_stats.py
Statistical summary script for VIROnator cohort master viral reports.

Reads a consolidated master viral report TSV (or individual report TSVs)
and generates:
1. cohort_stats_summary.tsv (9-column TSV)
2. cohort_stats_summary.md (Markdown summary snapshot)

Columns in cohort_stats_summary.tsv:
1. Phase
2. Project
3. Strategy
4. Total_Cohort_Samples
5. Positive_Samples
6. Prevalence_Pct
7. Total_Mapped_Reads
8. Top_Member_Name
9. Unique_Viruses_Detected
"""

import argparse
import collections
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Generate cohort viral statistical summary.")
    parser.add_argument("--input-report", required=True, help="Path to master viral report TSV file")
    parser.add_argument("--out-dir", default=".", help="Output directory for stats summary files (default: .)")
    return parser.parse_args()


def generate_stats(input_report_path, out_dir):
    if not os.path.exists(input_report_path):
        print(f"[ERROR] Input master report file not found: {input_report_path}")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    tsv_out_path = os.path.join(out_dir, "cohort_stats_summary.tsv")
    md_out_path = os.path.join(out_dir, "cohort_stats_summary.md")

    # Data structure: grouped by (phase, project, strategy)
    # Stores: { group_key: { "samples": set(), "pos_samples": set(), "total_reads": int, "virus_counts": Counter(), "unique_viruses": set() } }
    groups = collections.defaultdict(lambda: {
        "samples": set(),
        "pos_samples": set(),
        "total_reads": 0,
        "virus_counts": collections.Counter(),
        "unique_viruses": set()
    })

    standard_14_header = [
        "sample_id", "virus_accession", "virus_length", "virus_mapped_reads",
        "normalized_coverage", "physical_coverage", "human_genome_size",
        "sample_read_depth", "viral_copy_number", "virus_name_sanitized",
        "specimen", "phase", "project", "source_file"
    ]

    with open(input_report_path, "r") as f:
        header = None
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            parts = line_str.split("\t")
            
            if header is None:
                first_elem = parts[0].lower().strip()
                if "sample" in first_elem or "virus" in first_elem or "phase" in first_elem:
                    header = [p.lower().strip() for p in parts]
                    continue
                else:
                    # Line 1 is a data line without header: use standard 14-column header
                    header = standard_14_header[:len(parts)]

            row = dict(zip(header, parts))
            sample_id = row.get("sample_id", parts[0] if len(parts) > 0 else "").strip()
            virus_acc = row.get("virus_accession", parts[1] if len(parts) > 1 else "").strip()
            mapped_reads_str = row.get("virus_mapped_reads", parts[3] if len(parts) > 3 else "0").strip()
            virus_name = row.get("virus_name_sanitized", parts[9] if len(parts) > 9 else "").strip()
            specimen = row.get("specimen", parts[10] if len(parts) > 10 else "Unknown").strip()
            phase = row.get("phase", parts[11] if len(parts) > 11 else "unknown").strip()
            project = row.get("project", parts[12] if len(parts) > 12 else "base").strip()
            strategy = row.get("source_file", parts[13] if len(parts) > 13 else "unknown").strip()

            if not project:
                project = "base"

            try:
                mapped_reads = int(mapped_reads_str)
            except ValueError:
                mapped_reads = 0

            group_key = (phase, project, strategy)
            g_data = groups[group_key]
            
            g_data["samples"].add(sample_id)

            if virus_acc and virus_acc.lower() != "none" and mapped_reads > 0:
                g_data["pos_samples"].add(sample_id)
                g_data["total_reads"] += mapped_reads
                display_name = virus_name if virus_name and virus_name.lower() != "none" else virus_acc
                g_data["virus_counts"][display_name] += mapped_reads
                g_data["unique_viruses"].add(virus_acc)

    header_cols = [
        "Phase",
        "Project",
        "Strategy",
        "Total_Cohort_Samples",
        "Positive_Samples",
        "Prevalence_Pct",
        "Total_Mapped_Reads",
        "Top_Member_Name",
        "Unique_Viruses_Detected"
    ]

    summary_rows = []
    
    # Sort groups by Phase, Project, Strategy
    sorted_group_keys = sorted(groups.keys(), key=lambda x: (x[0], x[1], x[2]))

    for group_key in sorted_group_keys:
        phase, project, strategy = group_key
        g_data = groups[group_key]

        total_samples = len(g_data["samples"])
        pos_samples = len(g_data["pos_samples"])
        prev_pct = (pos_samples / float(total_samples) * 100.0) if total_samples > 0 else 0.0
        total_reads = g_data["total_reads"]
        unique_viruses_cnt = len(g_data["unique_viruses"])

        if g_data["virus_counts"]:
            top_member = g_data["virus_counts"].most_common(1)[0][0]
        else:
            top_member = "None"

        row = [
            phase,
            project,
            strategy,
            str(total_samples),
            str(pos_samples),
            f"{prev_pct:.2f}",
            str(total_reads),
            top_member,
            str(unique_viruses_cnt)
        ]
        summary_rows.append(row)

    # Write TSV summary
    with open(tsv_out_path, "w") as f:
        f.write("\t".join(header_cols) + "\n")
        for r in summary_rows:
            f.write("\t".join(r) + "\n")

    # Write Markdown summary
    with open(md_out_path, "w") as f:
        f.write("# VIROnator Cohort Statistical Summary Snapshot\n\n")
        f.write(f"**Input Master Report:** `{input_report_path}`  \n")
        f.write(f"**Total Sub-group Evaluation Combinations:** `{len(summary_rows)}`  \n\n")
        f.write("| " + " | ".join(header_cols) + " |\n")
        f.write("| " + " | ".join(["---"] * len(header_cols)) + " |\n")
        for r in summary_rows:
            f.write("| " + " | ".join(r) + " |\n")

    print(f"[SUCCESS] Generated 9-column stats summary TSV: {tsv_out_path}")
    print(f"[SUCCESS] Generated Markdown summary snapshot:   {md_out_path}")


def main():
    args = parse_args()
    generate_stats(args.input_report, args.out_dir)


if __name__ == "__main__":
    main()
