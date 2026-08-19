#!/usr/bin/env python3
"""
generate_distributions.py
Generates per-virus read distribution TIFF plots (plots/) and virus statistics TSV (stats/)
from consolidated master viral report TSV.

Outputs:
1. plots/dist_<PHASE>_<PROJECT>_<STRATEGY>_<VIRUS_ID>_<VIRUS_NAME>.tif (300 DPI TIFF)
2. stats/virus_stats_summary.tsv
"""

import argparse
import collections
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Generate per-virus read distribution plots and statistics.")
    parser.add_argument("--input-report", required=True, help="Path to master viral report TSV file")
    parser.add_argument("--out-dir", default=".", help="Base VIROnator directory (default: .)")
    parser.add_argument("--plots-dir", default=None, help="Specific output directory for TIFF plots")
    parser.add_argument("--stats-dir", default=None, help="Specific output directory for virus stats TSV")
    parser.add_argument("--target-phase", default=None, help="Target phase to filter (e.g. phase1 or 1)")
    parser.add_argument("--target-project", default=None, help="Target project to filter (e.g. base)")
    parser.add_argument("--strategies", nargs="*", default=[
        "exogeneSR_viral_clean_filtered.sorted.flags.cram",
        "exogeneSR_viral_raw_filtered.sorted.flags.cram"
    ], help="Target CRAM strategy file names to evaluate")
    return parser.parse_args()


def sanitize_filename(name):
    """Sanitizes strings for safe filename usage."""
    import re
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", name).strip("_")


def get_short_strategy(fname):
    """Maps long CRAM filename to concise strategy label for titles and filenames."""
    if "clean_filtered.sorted.flags" in fname:
        return "clean_flags"
    elif "raw_filtered.sorted.flags" in fname:
        return "raw_flags"
    elif "clean_filtered" in fname:
        return "clean_filtered"
    elif "raw_filtered" in fname:
        return "raw_filtered"
    elif "clean" in fname:
        return "clean_unfiltered"
    elif "raw" in fname:
        return "raw_unfiltered"
    return sanitize_filename(fname)


def main():
    args = parse_args()

    if not os.path.exists(args.input_report):
        print(f"[ERROR] Master report file not found: {args.input_report}")
        sys.exit(1)

    plots_dir = args.plots_dir if args.plots_dir else os.path.join(args.out_dir, "plots")
    stats_dir = args.stats_dir if args.stats_dir else os.path.join(args.out_dir, "stats")
    target_strategies = set(args.strategies)

    # If writing to /mnt/disks/staff mount point, use a safe local buffer folder first
    use_local_buffer = plots_dir.startswith("/mnt/disks") or stats_dir.startswith("/mnt/disks")
    if use_local_buffer:
        actual_plots_dir = plots_dir
        actual_stats_dir = stats_dir
        plots_dir = os.path.abspath("plots_temp")
        stats_dir = os.path.abspath("stats_temp")

    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(stats_dir, exist_ok=True)

    # Phase & Project target normalization
    target_phase = str(args.target_phase).strip().lower() if args.target_phase else None
    if target_phase:
        if not target_phase.startswith("phase") and target_phase.isdigit():
            target_phase_alt = f"phase{target_phase}"
        else:
            target_phase_alt = target_phase.replace("phase", "")
    else:
        target_phase_alt = None

    target_project = str(args.target_project).strip().lower() if args.target_project else "base"
    if not target_project:
        target_project = "base"

    # Global tracking per dataset (phase, project, strategy)
    # total_samples_map: { (phase, project, strategy): set_of_all_samples }
    # viral_pos_samples_map: { (phase, project, strategy): set_of_virus_pos_samples }
    # virus_data: { (phase, project, strategy, virus_acc): { "name": str, "samples": dict(sample_id -> mapped_reads) } }

    total_samples_map = collections.defaultdict(set)
    viral_pos_samples_map = collections.defaultdict(set)
    virus_data = collections.defaultdict(lambda: {
        "name": "",
        "samples": {}
    })

    standard_14_header = [
        "sample_id", "virus_accession", "virus_length", "virus_mapped_reads",
        "normalized_coverage", "physical_coverage", "human_genome_size",
        "sample_read_depth", "viral_copy_number", "virus_name_sanitized",
        "specimen", "phase", "project", "source_file"
    ]

    with open(args.input_report, "r") as f:
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
                    header = standard_14_header[:len(parts)]

            row = dict(zip(header, parts))
            sample_id = row.get("sample_id", parts[0] if len(parts) > 0 else "").strip()
            virus_acc = row.get("virus_accession", parts[1] if len(parts) > 1 else "").strip()
            mapped_reads_str = row.get("virus_mapped_reads", parts[3] if len(parts) > 3 else "0").strip()
            virus_name = row.get("virus_name_sanitized", parts[9] if len(parts) > 9 else "").strip()
            phase = row.get("phase", parts[11] if len(parts) > 11 else "unknown").strip()
            project = row.get("project", parts[12] if len(parts) > 12 else "base").strip()
            strategy = row.get("source_file", parts[13] if len(parts) > 13 else "unknown").strip()

            if not project:
                project = "base"

            # Filter by requested target strategies if specified
            if target_strategies and strategy not in target_strategies:
                continue

            # Filter by requested target phase and project if specified
            row_phase_lower = phase.lower().strip()
            row_proj_lower = project.lower().strip()
            if not row_proj_lower:
                row_proj_lower = "base"

            if target_phase and not (row_phase_lower == target_phase or row_phase_lower == target_phase_alt):
                continue
            if target_project and row_proj_lower != target_project:
                continue

            try:
                mapped_reads = int(mapped_reads_str)
            except ValueError:
                mapped_reads = 0

            ds_key = (phase, project, strategy)
            total_samples_map[ds_key].add(sample_id)

            if virus_acc and virus_acc.lower() != "none" and mapped_reads > 0:
                viral_pos_samples_map[ds_key].add(sample_id)
                
                v_key = (phase, project, strategy, virus_acc)
                display_name = virus_name if virus_name and virus_name.lower() != "none" else virus_acc
                virus_data[v_key]["name"] = display_name
                virus_data[v_key]["samples"][sample_id] = mapped_reads

    # --------------------------------------------------------------------------
    # 1. Generate virus_stats_summary_<PHASE>_<PROJECT>.tsv in stats_dir
    # --------------------------------------------------------------------------
    phase_tag = f"phase{args.target_phase}" if args.target_phase and not str(args.target_phase).startswith("phase") else (str(args.target_phase) if args.target_phase else "all")
    proj_tag = str(args.target_project) if args.target_project else "base"
    stats_tsv_name = f"virus_stats_summary_{phase_tag}_{proj_tag}.tsv"
    stats_tsv_path = os.path.join(stats_dir, stats_tsv_name)
    stats_header = [
        "Phase",
        "Project",
        "Strategy",
        "Virus_Accession",
        "Virus_Name_Sanitized",
        "Positive_Samples_Count",
        "Total_Viral_Positive_Samples",
        "Total_Cohort_Samples",
        "Pct_Of_Viral_Positive_Samples",
        "Pct_Of_Total_Cohort_Samples",
        "Total_Reads_Assigned"
    ]

    stats_rows = []
    sorted_v_keys = sorted(virus_data.keys(), key=lambda x: (x[0], x[1], x[2], x[3]))

    for v_key in sorted_v_keys:
        phase, project, strategy, virus_acc = v_key
        ds_key = (phase, project, strategy)
        
        v_info = virus_data[v_key]
        sample_dict = v_info["samples"]
        
        pos_count = len(sample_dict)
        total_viral_pos = len(viral_pos_samples_map[ds_key])
        total_cohort = len(total_samples_map[ds_key])
        
        pct_viral_pos = (pos_count / float(total_viral_pos) * 100.0) if total_viral_pos > 0 else 0.0
        pct_cohort = (pos_count / float(total_cohort) * 100.0) if total_cohort > 0 else 0.0
        total_reads = sum(sample_dict.values())
        
        short_strat = get_short_strategy(strategy)

        row = [
            phase,
            project,
            short_strat,
            virus_acc,
            v_info["name"],
            str(pos_count),
            str(total_viral_pos),
            str(total_cohort),
            f"{pct_viral_pos:.2f}",
            f"{pct_cohort:.2f}",
            str(total_reads)
        ]
        stats_rows.append(row)

    with open(stats_tsv_path, "w") as f:
        f.write("\t".join(stats_header) + "\n")
        for r in stats_rows:
            f.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Written virus stats summary TSV: {stats_tsv_path}")

    # --------------------------------------------------------------------------
    # 2. Generate 2-Panel Publication Quality TIFF Plots in plots/
    # --------------------------------------------------------------------------
    plot_count = 0
    color_map = {
        "clean_filtered": "#80ed99",
        "clean_flags": "#57cc99",
        "raw_filtered": "#38a3a5",
        "raw_flags": "#22577a",
        "clean_unfiltered": "#80ed99",
        "raw_unfiltered": "#38a3a5"
    }

    from matplotlib.ticker import MaxNLocator

    total_plot_keys = len(sorted_v_keys)
    print(f"[INFO] Rendering {total_plot_keys} distribution plots...")

    for idx, v_key in enumerate(sorted_v_keys, 1):
        phase, project, strategy, virus_acc = v_key
        ds_key = (phase, project, strategy)
        
        v_info = virus_data[v_key]
        read_counts = sorted(v_info["samples"].values())
        if not read_counts:
            continue

        pos_count = len(read_counts)
        total_viral_pos = len(viral_pos_samples_map[ds_key])
        total_cohort = len(total_samples_map[ds_key])

        pct_viral_pos = (pos_count / float(total_viral_pos) * 100.0) if total_viral_pos > 0 else 0.0
        pct_cohort = (pos_count / float(total_cohort) * 100.0) if total_cohort > 0 else 0.0

        short_strat = get_short_strategy(strategy)
        v_name_sanitized = sanitize_filename(v_info["name"])

        out_tif = os.path.join(plots_dir, f"dist_{phase}_{project}_{short_strat}_{virus_acc}_{v_name_sanitized}.tiff")
        print(f"  [{idx}/{total_plot_keys}] Generating plot: {os.path.basename(out_tif)} ...", flush=True)

        # Publication Figure Setup (12 x 5 inches, 300 DPI)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5.2), dpi=300)
        main_color = color_map.get(short_strat, "#57cc99")

        max_r = max(read_counts)
        median_r = float(np.median(read_counts))
        q25_r = float(np.percentile(read_counts, 25))
        q75_r = float(np.percentile(read_counts, 75))

        from matplotlib.ticker import MaxNLocator, FormatStrFormatter

        # ----------------------------------------------------------------------
        # Panel A: Read Count Frequency Distribution (Strict Integer X-Axis & Spaced Bars)
        # ----------------------------------------------------------------------
        unique_vals, val_counts = np.unique(read_counts, return_counts=True)
        bar_width = 0.35  # Thin bars with 65% clear spacing between bars
        ax1.bar(unique_vals, val_counts, width=bar_width, color=main_color, edgecolor='#ffffff', linewidth=0.8, alpha=0.9, align='center')
        ax1.set_xlabel("Mapped Read Count", fontsize=11, fontweight='bold', labelpad=8, color='black')
        
        # Enforce pure integer ticks (no decimals)
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        if max_r <= 25:
            ax1.set_xticks(np.arange(min(unique_vals), max(unique_vals) + 1, dtype=int))
        else:
            ax1.set_xticks(unique_vals)

        ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax1.set_ylabel("Number of Samples", fontsize=11, fontweight='bold', labelpad=8, color='black')
        ax1.set_title("A. Read Count Frequency Distribution", fontsize=12, fontweight='bold', pad=12, color='black')
        
        # Despine top and right axes
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        # ----------------------------------------------------------------------
        # Panel B: Ordered Read Counts (Solid-filled dots, no median line, no legend)
        # ----------------------------------------------------------------------
        ranks = np.arange(1, pos_count + 1)
        ax2.plot(ranks, read_counts, marker='o', color=main_color, linewidth=2, markersize=6, markerfacecolor=main_color, markeredgecolor=main_color)
        ax2.fill_between(ranks, read_counts, color=main_color, alpha=0.15)
        
        # Sleek Summary Callout Badge
        badge_text = f"Max Reads: {max_r}\nMedian Reads: {median_r:.1f}\nIQR: [{q25_r:.1f} - {q75_r:.1f}]"
        ax2.text(0.04, 0.94, badge_text, transform=ax2.transAxes, fontsize=9.5, verticalalignment='top',
                 bbox=dict(boxstyle='round,pad=0.5', facecolor='#f7fafc', edgecolor='#cbd5e0', alpha=0.92))

        ax2.xaxis.set_major_locator(MaxNLocator(integer=True, steps=[1, 2, 5, 10]))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        if pos_count <= 25:
            ax2.set_xticks(ranks)
        else:
            ax2.set_xticks(np.linspace(1, pos_count, min(10, pos_count), dtype=int))
        ax2.set_xlabel("Sample Occurrence Rank (Increasing Order)", fontsize=11, fontweight='bold', labelpad=8, color='black')
        ax2.set_ylabel("Mapped Read Count", fontsize=11, fontweight='bold', labelpad=8, color='black')
        ax2.set_title("B. Ordered Read Counts Across Samples", fontsize=12, fontweight='bold', pad=12, color='black')

        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        # Main Title (Running across both panels, no angled brackets)
        main_title = f"Phase: {phase} | Project: {project} | Strategy: {short_strat}\n{v_info['name']} ({virus_acc}) N = {pos_count} ({pct_viral_pos:.2f}%* / {pct_cohort:.2f}%**)"
        fig.suptitle(main_title, fontsize=13, fontweight='bold', y=0.98, color='black')

        # Subtitle Footnote Legend (Clean black text)
        footnote = f"* % of total virus-positive samples in dataset ({pos_count} / {total_viral_pos})     ** % of total samples in dataset ({pos_count} / {total_cohort})"
        fig.text(0.5, 0.012, footnote, ha='center', fontsize=9, style='italic', color='black')

        plt.tight_layout(rect=[0, 0.06, 1, 0.91])

        plt.savefig(out_tif, format='tiff', dpi=300)
        plt.close(fig)
        plot_count += 1

    print(f"[SUCCESS] Generated {plot_count} high-resolution TIFF distribution plots in {plots_dir}")

    if use_local_buffer:
        import shutil, subprocess
        print(f"[INFO] Syncing plots and stats to target destination: {actual_plots_dir} and {actual_stats_dir} ...")
        # Try local copy first, fallback to gsutil cp
        try:
            subprocess.run(f"mkdir -p '{actual_plots_dir}' '{actual_stats_dir}' 2>/dev/null", shell=True)
            for f in os.listdir(stats_dir):
                shutil.copy(os.path.join(stats_dir, f), os.path.join(actual_stats_dir, f))
            for f in os.listdir(plots_dir):
                shutil.copy(os.path.join(plots_dir, f), os.path.join(actual_plots_dir, f))
            print("[SUCCESS] Successfully copied all output files to staff mount point.")
        except Exception:
            bucket_name = "ml-phi-staff-m277455-p-rsa-us-central1-p-a3d4"
            subprocess.run(f"gsutil -q -m cp -r {stats_dir}/* gs://{bucket_name}/SSC_hg38_stats/ 2>/dev/null", shell=True)
            subprocess.run(f"gsutil -q -m cp -r {plots_dir}/* gs://{bucket_name}/SSC_hg38_plots/ 2>/dev/null", shell=True)
            print("[SUCCESS] Successfully synced all output files to GCS bucket via gsutil.")


if __name__ == "__main__":
    main()
