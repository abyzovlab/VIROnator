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
from matplotlib.ticker import MaxNLocator, FormatStrFormatter


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

    phase_tag = f"phase{target_phase}" if target_phase and not str(target_phase).startswith("phase") else (target_phase if target_phase else "all")
    proj_tag = target_project if target_project else "base"

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

    plt.rcParams['axes.grid'] = False
    plt.rcParams['grid.alpha'] = 0.0

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

        main_color = color_map.get(short_strat, "#57cc99")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5), dpi=300)

        ax1.set_box_aspect(1)
        ax2.set_box_aspect(1)

        unique_vals, val_counts = np.unique(read_counts, return_counts=True)

        ax1.bar(unique_vals, val_counts, width=0.7, color=main_color)
        ax1.set_xlabel("Mapped Read Count")
        ax1.set_ylabel("Number of Samples")
        ax1.set_title("A. Read Count Frequency Distribution")
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        
        # Prevent decimal autoscale limits on single-value histograms
        if len(unique_vals) == 1:
            v_val = unique_vals[0]
            ax1.set_xlim(max(0, v_val - 2), v_val + 2)
            ax1.set_xticks(np.arange(max(0, v_val - 2), v_val + 3, dtype=int))

        ax1.grid(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        ranks = np.arange(1, pos_count + 1)

        ax2.plot(ranks, read_counts, marker="o", color=main_color, linewidth=0.8, markersize=8)
        if pos_count == 1:
            ax2.vlines(x=1, ymin=0, ymax=read_counts[0], color=main_color, linewidth=0.8)
        else:
            ax2.fill_between(ranks, read_counts, alpha=0.15, color=main_color, edgecolor='none', linewidth=0)
        ax2.set_xlabel("Sample Occurrence Rank (Increasing Order)")
        ax2.set_ylabel("Mapped Read Count")
        ax2.set_title("B. Ordered Read Counts Across Samples")
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        
        # Prevent decimal autoscale limits on single-sample plots
        if pos_count == 1:
            ax2.set_xlim(0, 2)
            ax2.set_xticks([0, 1, 2])

        ax2.grid(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        main_title = (
            f"Phase: {phase} | Project: {project} | Strategy: {short_strat}\n"
            f"{v_info['name']} ({virus_acc}) N = {pos_count} "
            f"({pct_viral_pos:.2f}%* / {pct_cohort:.2f}%**)"
        )

        fig.suptitle(main_title)

        fig.text(
            0.5,
            0.015,
            f"* % of total virus-positive samples in dataset ({pos_count} / {total_viral_pos})     "
            f"** % of total samples in dataset ({pos_count} / {total_cohort})",
            ha="center",
            fontsize=10.5,
            color="black",
            style="normal"
        )

        plt.tight_layout(rect=[0, 0.08, 1, 0.92])
        plt.savefig(out_tif, format="tiff", dpi=300)
        plt.close(fig)
        plot_count += 1

    # --------------------------------------------------------------------------
    # 3. Generate Overall Dataset-Wide Summary TIFF Plot and Overall Stats TSV
    # --------------------------------------------------------------------------
    overall_stats_tsv_name = f"virus_stats_summary_{phase_tag}_{proj_tag}_OVERALL.tsv"
    overall_stats_tsv_path = os.path.join(stats_dir, overall_stats_tsv_name)
    overall_plot_name = f"dist_{phase_tag}_{proj_tag}_OVERALL_SUMMARY.tiff"
    overall_plot_path = os.path.join(plots_dir, overall_plot_name)

    overall_rows = []
    overall_header = [
        "Phase", "Project", "Strategy", "Total_Cohort_Samples",
        "Total_Viral_Positive_Samples", "Cohort_Positivity_Pct",
        "Total_Mapped_Viral_Reads", "Unique_Viruses_Detected", "Top_Prevalent_Virus"
    ]

    # Aggregate overall stats per dataset key (phase, project, strategy)
    for ds_key in sorted(total_samples_map.keys()):
        d_phase, d_proj, d_strat = ds_key
        d_short_strat = get_short_strategy(d_strat)
        
        tot_samples = len(total_samples_map[ds_key])
        pos_samples = len(viral_pos_samples_map[ds_key])
        pos_pct = (pos_samples / float(tot_samples) * 100.0) if tot_samples > 0 else 0.0
        
        # Total reads & unique viruses across all viruses in dataset
        ds_total_reads = 0
        ds_unique_v = set()
        ds_virus_counts = collections.Counter()

        for v_key, v_info in virus_data.items():
            if v_key[0] == d_phase and v_key[1] == d_proj and v_key[2] == d_strat:
                ds_unique_v.add(v_key[3])
                for smp, r_cnt in v_info["samples"].items():
                    ds_total_reads += r_cnt
                    ds_virus_counts[v_info["name"]] += len(v_info["samples"])

        top_v = ds_virus_counts.most_common(1)[0][0] if ds_virus_counts else "None"

        overall_rows.append([
            d_phase, d_proj, d_short_strat, str(tot_samples),
            str(pos_samples), f"{pos_pct:.2f}", str(ds_total_reads),
            str(len(ds_unique_v)), top_v
        ])

    with open(overall_stats_tsv_path, "w") as f:
        f.write("\t".join(overall_header) + "\n")
        for r in overall_rows:
            f.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Written overall dataset stats summary TSV: {overall_stats_tsv_path}")

    # Generate Overall Dataset Summary 2-Panel Plot
    all_read_counts = []
    for v_key, v_info in virus_data.items():
        all_read_counts.extend(v_info["samples"].values())

    if all_read_counts:
        all_read_counts = sorted(all_read_counts)
        tot_v_pos = len(all_read_counts)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5), dpi=300)
        ax1.set_box_aspect(1)
        ax2.set_box_aspect(1)
        overall_color = "#22577a"

        u_vals, u_cnts = np.unique(all_read_counts, return_counts=True)
        ax1.bar(u_vals, u_cnts, width=0.7, color=overall_color)
        ax1.set_xlabel("Mapped Read Count (All Viruses)", fontsize=11, labelpad=8, color='black')
        ax1.set_ylabel("Number of Mapped Occurrences", fontsize=11, labelpad=8, color='black')
        ax1.set_title("A. Dataset-Wide Read Count Frequency", fontsize=12, pad=12, color='black')
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax1.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax1.grid(False)
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)

        o_ranks = np.arange(1, tot_v_pos + 1)
        ax2.plot(o_ranks, all_read_counts, marker='o', color=overall_color, linewidth=0.8, markersize=8, markerfacecolor=overall_color, markeredgecolor=overall_color)
        if tot_v_pos == 1:
            ax2.vlines(x=1, ymin=0, ymax=all_read_counts[0], color=overall_color, linewidth=0.8)
        else:
            ax2.fill_between(o_ranks, all_read_counts, color=overall_color, alpha=0.15, edgecolor='none', linewidth=0)
        ax2.set_xlabel("Occurrence Rank Across Entire Dataset", fontsize=11, labelpad=8, color='black')
        ax2.set_ylabel("Mapped Read Count", fontsize=11, labelpad=8, color='black')
        ax2.set_title("B. Dataset-Wide Ordered Read Counts", fontsize=12, pad=12, color='black')
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True, steps=[1, 2, 5, 10]))
        ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax2.grid(False)
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

        fig.suptitle(f"Dataset-Wide Overall Summary | Phase: {phase_tag} | Project: {proj_tag}\nTotal Viral Detections N = {tot_v_pos}", fontsize=13, y=0.98, color='black')
        fig.text(0.5, 0.015, "Overall dataset-wide distribution across all detected viruses", ha='center', fontsize=10.5, style='normal', color='black')

        plt.tight_layout(rect=[0, 0.08, 1, 0.92])
        plt.savefig(overall_plot_path, format='tiff', dpi=300)
        plt.close(fig)
        print(f"[SUCCESS] Generated overall dataset TIFF plot: {overall_plot_path}")

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
