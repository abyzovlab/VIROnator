import argparse
import collections
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator, FormatStrFormatter


def parse_args():
    parser = argparse.ArgumentParser(description="Generate cohort viral statistical summary and plots.")
    parser.add_argument("--input-report", required=True, help="Path to master viral report TSV file")
    parser.add_argument("--out-dir", default=".", help="Output directory for stats summary files (default: .)")
    parser.add_argument("--dataset", default="cohort", help="Dataset identifier (e.g. SSC, MCBiobank)")
    parser.add_argument("--genome-build", default="hg38", help="Genome build identifier (e.g. hg38)")
    parser.add_argument("--target-phase", default=None, help="Target phase to filter")
    parser.add_argument("--target-project", default=None, help="Target project to filter")
    parser.add_argument("--strategies", nargs="*", default=[
        "exogeneSR_viral_clean_filtered.sorted.flags.cram"
    ], help="Target CRAM strategy file names")
    parser.add_argument("--cohort-scope", default="combined_all", choices=["target_only", "combined_all", "both"],
                        help="Cohort scope setting")
    parser.add_argument("--reads-panel-a-loglog", default="on", choices=["on", "off"],
                        help="Reads distribution Panel A scale: 'on' for Log-Log, 'off' for linear")
    parser.add_argument("--reads-panel-b-log-y", default="on", choices=["on", "off"],
                        help="Reads distribution Panel B scale: 'on' for Log Y, 'off' for linear")
    parser.add_argument("--copy-number-panel-a-loglog", default="on", choices=["on", "off"],
                        help="Copy number distribution Panel A scale: 'on' for Log-Log, 'off' for linear")
    parser.add_argument("--copy-number-panel-b-log-y", default="on", choices=["on", "off"],
                        help="Copy number distribution Panel B scale: 'on' for Log Y, 'off' for linear")
    parser.add_argument("--log-scale-read-cutoff", type=int, default=30,
                        help="Read count threshold to trigger Log-Scale transformations")
    parser.add_argument("--log-scale-copy-number-cutoff", type=float, default=0.1,
                        help="Copy number threshold to trigger Log-Scale transformations")
    parser.add_argument("--prelim-prevalence-cutoff-pct", type=float, default=5.0,
                        help="Preliminary prevalence cutoff %")
    parser.add_argument("--prelim-mean-read-cutoff", type=float, default=6.0,
                        help="Preliminary mean read count cutoff")
    parser.add_argument("--heatmap-read-counts", default="off", choices=["on", "off"],
                        help="Generate read counts heatmap in stats directory")
    parser.add_argument("--heatmap-copy-number", default="off", choices=["on", "off"],
                        help="Generate copy number heatmap in stats directory")
    parser.add_argument("--target-heatmap-strategy", default="clean_flags",
                        help="Strategy filter for heatmap generation")
    return parser.parse_args()


def classify_virus(pct_cohort, mean_reads, prev_cutoff=5.0, mean_cutoff=6.0):
    if mean_reads <= mean_cutoff:
        if pct_cohort < prev_cutoff:
            return "sporadic_noise"
        else:
            return "systematic_noise"
    else:
        if pct_cohort >= prev_cutoff:
            return "virome"
        else:
            return "infection"


def get_short_strategy(fname):
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
    import re
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", fname).strip("_")


def generate_stats(args):
    input_report_path = args.input_report
    out_dir = args.out_dir
    dataset = args.dataset or "cohort"
    genome_build = args.genome_build or "hg38"

    if not os.path.exists(input_report_path):
        print(f"[ERROR] Input master report file not found: {input_report_path}")
        sys.exit(1)

    os.makedirs(out_dir, exist_ok=True)
    prefix = f"{dataset}_{genome_build}"
    tsv_out_path = os.path.join(out_dir, f"{prefix}_stats_summary.tsv")
    md_out_path = os.path.join(out_dir, f"{prefix}_stats_summary.md")

    target_strategies = set(args.strategies) if args.strategies else None

    # Normalization of phase & project
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

    cohort_scope = args.cohort_scope.strip().lower() if args.cohort_scope else "combined_all"

    # Data structures for master summary and per-virus stats
    total_samples_map = collections.defaultdict(set)
    viral_pos_samples_map = collections.defaultdict(set)
    virus_data = collections.defaultdict(lambda: {
        "name": "",
        "samples": dict(), # sample_id -> mapped_reads
        "copy_numbers": dict() # sample_id -> copy_number
    })

    # Record-level master data for reads & copy_number overall summary plots
    dataset_records = collections.defaultdict(lambda: {"reads": [], "copy_numbers": []})

    standard_14_header = [
        "sample_id", "virus_accession", "virus_length", "virus_mapped_reads",
        "normalized_coverage", "physical_coverage", "human_genome_size",
        "sample_read_depth", "viral_copy_number", "virus_name_sanitized",
        "specimen", "phase", "project", "source_file"
    ]

    # Fast vector reading using pandas
    try:
        df_raw = pd.read_csv(input_report_path, sep="\t", dtype=str)
    except Exception as e:
        print(f"[ERROR] Failed to read report TSV with pandas: {e}")
        sys.exit(1)

    df_raw.columns = [c.lower().strip() for c in df_raw.columns]
    
    # Required columns fallback
    col_map = {
        'sample_id': df_raw.columns[0] if len(df_raw.columns) > 0 else 'sample_id',
        'virus_accession': df_raw.columns[1] if len(df_raw.columns) > 1 else 'virus_accession',
        'virus_mapped_reads': df_raw.columns[3] if len(df_raw.columns) > 3 else 'virus_mapped_reads',
        'viral_copy_number': df_raw.columns[8] if len(df_raw.columns) > 8 else 'viral_copy_number',
        'virus_name_sanitized': df_raw.columns[9] if len(df_raw.columns) > 9 else 'virus_name_sanitized',
        'phase': df_raw.columns[11] if len(df_raw.columns) > 11 else 'phase',
        'project': df_raw.columns[12] if len(df_raw.columns) > 12 else 'project',
        'source_file': df_raw.columns[13] if len(df_raw.columns) > 13 else 'source_file'
    }

    df = pd.DataFrame()
    df['sample_id'] = df_raw[col_map['sample_id']].fillna('').astype(str).str.strip()
    df['virus_accession'] = df_raw[col_map['virus_accession']].fillna('').astype(str).str.strip()
    df['virus_mapped_reads'] = pd.to_numeric(df_raw[col_map['virus_mapped_reads']], errors='coerce').fillna(0).astype(int)
    df['viral_copy_number'] = pd.to_numeric(df_raw[col_map['viral_copy_number']], errors='coerce').fillna(0.0).astype(float)
    df['virus_name_sanitized'] = df_raw[col_map['virus_name_sanitized']].fillna('').astype(str).str.strip()
    df['phase'] = df_raw[col_map['phase']].fillna('unknown').astype(str).str.strip()
    df['project'] = df_raw[col_map['project']].fillna('base').astype(str).str.strip()
    df['project'] = df['project'].replace({'': 'base', 'none': 'base', '0': 'base'})
    df['source_file'] = df_raw[col_map['source_file']].fillna('unknown').astype(str).str.strip()

    if target_strategies:
        df = df[df['source_file'].isin(target_strategies)].copy()

    # Normalize phase matching
    df['phase_clean'] = df['phase'].str.lower().str.strip()
    df['project_clean'] = df['project'].str.lower().str.strip()

    # Pre-build data structures expected downstream
    total_samples_map = collections.defaultdict(set)
    viral_pos_samples_map = collections.defaultdict(set)
    virus_data = collections.defaultdict(lambda: {
        "name": "",
        "samples": dict(),
        "copy_numbers": dict()
    })
    dataset_records = collections.defaultdict(lambda: {"reads": [], "copy_numbers": []})

    # Vector iteration using itertuples for speed
    for row in df.itertuples(index=False):
        sample_id = row.sample_id
        virus_acc = row.virus_accession
        mapped_reads = row.virus_mapped_reads
        copy_num = row.viral_copy_number
        virus_name = row.virus_name_sanitized
        raw_phase = row.phase
        raw_proj_clean = row.project_clean
        strategy = row.source_file

        raw_phase_clean = row.phase_clean

        matches_target = True
        if target_phase and not (raw_phase_clean == target_phase or raw_phase_clean == target_phase_alt):
            matches_target = False
        if target_project and raw_proj_clean != target_project:
            matches_target = False

        keys_to_add = []
        if cohort_scope in ["target_only", "both"] and matches_target:
            keys_to_add.append((raw_phase, raw_proj_clean))
        if cohort_scope in ["combined_all", "both"]:
            keys_to_add.append(("all_cohorts", "combined"))

        for p_val, prj_val in keys_to_add:
            ds_key = (p_val, prj_val, strategy)
            total_samples_map[ds_key].add(sample_id)

            if virus_acc and virus_acc.lower() != "none" and mapped_reads > 0:
                viral_pos_samples_map[ds_key].add(sample_id)
                v_key = (p_val, prj_val, strategy, virus_acc)
                display_name = virus_name if virus_name and virus_name.lower() != "none" else virus_acc
                virus_data[v_key]["name"] = display_name
                virus_data[v_key]["samples"][sample_id] = mapped_reads
                virus_data[v_key]["copy_numbers"][sample_id] = copy_num

                ds_group_key = (p_val, prj_val)
                dataset_records[ds_group_key]["reads"].append(mapped_reads)
                dataset_records[ds_group_key]["copy_numbers"].append(copy_num)

    # 1. Standard 9-column Cohort Stats Summary TSV & MD
    header_cols = [
        "Phase", "Project", "Strategy", "Total_Cohort_Samples",
        "Positive_Samples", "Prevalence_Pct", "Total_Mapped_Reads",
        "Top_Member_Name", "Unique_Viruses_Detected"
    ]
    summary_rows = []
    for ds_key in sorted(total_samples_map.keys()):
        p_val, prj_val, strat = ds_key
        tot_s = len(total_samples_map[ds_key])
        pos_s = len(viral_pos_samples_map[ds_key])
        prev_p = (pos_s / float(tot_s) * 100.0) if tot_s > 0 else 0.0
        
        relevant_vkeys = [v_k for v_k in virus_data.keys() if v_k[0] == p_val and v_k[1] == prj_val and v_k[2] == strat]
        tot_r = sum(sum(virus_data[v_k]["samples"].values()) for v_k in relevant_vkeys)
        u_vir = len({v_k[3] for v_k in relevant_vkeys})

        v_counts = collections.Counter()
        for v_k in relevant_vkeys:
            v_counts[virus_data[v_k]["name"]] += len(virus_data[v_k]["samples"])
        top_m = v_counts.most_common(1)[0][0] if v_counts else "None"

        short_s = get_short_strategy(strat)
        summary_rows.append([p_val, prj_val, short_s, str(tot_s), str(pos_s), f"{prev_p:.2f}", str(tot_r), top_m, str(u_vir)])

    with open(tsv_out_path, "w") as f:
        f.write("\t".join(header_cols) + "\n")
        for r in summary_rows:
            f.write("\t".join(r) + "\n")

    with open(md_out_path, "w") as f:
        f.write("# VIROnator Cohort Statistical Summary Snapshot\n\n")
        f.write(f"**Input Master Report:** `{input_report_path}`  \n\n")
        f.write("| " + " | ".join(header_cols) + " |\n")
        f.write("| " + " | ".join(["---"] * len(header_cols)) + " |\n")
        for r in summary_rows:
            f.write("| " + " | ".join(r) + " |\n")

    print(f"[SUCCESS] Written cohort stats summary TSV: {tsv_out_path}")
    print(f"[SUCCESS] Written cohort stats summary MD:  {md_out_path}")

    # Group virus data by (phase, project)
    dataset_groups = collections.defaultdict(list)
    for v_key in virus_data.keys():
        dataset_groups[(v_key[0], v_key[1])].append(v_key)

    for (cur_phase, cur_proj), group_v_keys in dataset_groups.items():
        is_phase_empty = not cur_phase or str(cur_phase).lower() in ["none", "0", "", "all_cohorts"]
        is_proj_empty = not cur_proj or str(cur_proj).lower() in ["none", "0", "", "base", "combined"]

        if is_phase_empty and is_proj_empty:
            p_file_tag = ""
            prj_file_tag = ""
            cur_phase_tag = "all"
            cur_proj_tag = "all"
            fn_prefix = f"{dataset}_{genome_build}_{strat_tag}"
        else:
            if is_phase_empty:
                p_file_tag = "all" if str(cur_phase).lower() == "all_cohorts" else "phase"
                cur_phase_tag = "all"
            elif str(cur_phase).startswith("phase"):
                p_file_tag = str(cur_phase)
                cur_phase_tag = str(cur_phase)
            else:
                p_file_tag = f"phase{cur_phase}"
                cur_phase_tag = f"phase{cur_phase}"

            if is_proj_empty:
                prj_file_tag = "combined" if str(cur_proj).lower() == "combined" else "base"
                cur_proj_tag = "base"
            else:
                prj_file_tag = str(cur_proj)
                cur_proj_tag = str(cur_proj)

            fn_prefix = f"{dataset}_{genome_build}_{p_file_tag}_{prj_file_tag}_{strat_tag}"

        # ----------------------------------------------------------------------
        # 1. Generate VIRUSES TSV & VIRUSES_classes.tsv
        # ----------------------------------------------------------------------
        viruses_tsv_name = f"{fn_prefix}_virus_stats_VIRUSES.tsv"
        viruses_tsv_path = os.path.join(out_dir, viruses_tsv_name)
        classes_tsv_name = f"{fn_prefix}_virus_stats_VIRUSES_classes.tsv"
        classes_tsv_path = os.path.join(out_dir, classes_tsv_name)

        viruses_header = [
            "Phase", "Project", "Strategy", "Virus_Accession", "Virus_Name_Sanitized",
            "Positive_Samples_Count", "Total_Viral_Positive_Samples", "Total_Cohort_Samples",
            "Pct_Of_Viral_Positive_Samples", "Pct_Of_Total_Cohort_Samples",
            "Total_Reads_Assigned", "Mean_Mapped_Reads_Per_Positive_Sample",
            "Preliminary_Classification"
        ]

        viruses_rows = []
        class_counts = collections.Counter({
            "systematic_noise": 0,
            "sporadic_noise": 0,
            "infection": 0,
            "virome": 0
        })

        sorted_group_v_keys = sorted(group_v_keys, key=lambda x: len(virus_data[x]["samples"]), reverse=True)

        for v_key in sorted_group_v_keys:
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
            mean_reads_per_pos = (total_reads / float(pos_count)) if pos_count > 0 else 0.0

            prelim_class = classify_virus(pct_cohort, mean_reads_per_pos, args.prelim_prevalence_cutoff_pct, args.prelim_mean_read_cutoff)
            class_counts[prelim_class] += 1
            short_strat = get_short_strategy(strategy)

            row = [
                phase, project, short_strat, virus_acc, v_info["name"],
                str(pos_count), str(total_viral_pos), str(total_cohort),
                f"{pct_viral_pos:.2f}", f"{pct_cohort:.2f}",
                str(total_reads), f"{mean_reads_per_pos:.2f}", prelim_class
            ]
            viruses_rows.append(row)

        with open(viruses_tsv_path, "w") as f:
            f.write("\t".join(viruses_header) + "\n")
            for r in viruses_rows:
                f.write("\t".join(r) + "\n")
        print(f"[SUCCESS] Written VIRUSES stats summary TSV: {viruses_tsv_path}")

        # ----------------------------------------------------------------------
        # 1b. VIRUSES_classes.tiff (Vertical bar plot of 4 classification categories ordered highest to lowest)
        # ----------------------------------------------------------------------
        classes_plot_name = f"{fn_prefix}_virus_stats_VIRUSES_classes.tiff"
        classes_plot_path = os.path.join(out_dir, classes_plot_name)

        # Sort classes descending by count
        sorted_classes = class_counts.most_common() # [(cat, count), ...]
        cls_names = [c[0] for c in sorted_classes]
        cls_counts = [c[1] for c in sorted_classes]

        fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
        bars = ax.bar(cls_names, cls_counts, width=0.88, color="#008fbf", edgecolor='none')
        ax.set_ylabel("Occurrence Count", fontsize=14, color='black', labelpad=8)
        ax.set_xlabel("Preliminary Classification Category", fontsize=14, color='black', labelpad=8)
        ax.tick_params(axis='x', labelsize=12, rotation=15)
        ax.tick_params(axis='y', labelsize=12)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_title(
            f"{dataset} Preliminary Classification Category Counts\n"
            f"Phase: {cur_phase_tag} | Project: {cur_proj_tag}",
            fontsize=16, pad=12, color='black'
        )

        max_c = max(cls_counts) if cls_counts else 1
        ax.set_ylim(0, max_c * 1.18)

        for bar in bars:
            h = bar.get_height()
            ax.annotate(f"{int(h)}", xy=(bar.get_x() + bar.get_width() / 2, h),
                        xytext=(0, 4), textcoords="offset points", ha='center', va='bottom', fontsize=12, color='black')

        ax.grid(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.savefig(classes_plot_path, format="tiff", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[SUCCESS] Generated VIRUSES classes TIFF plot: {classes_plot_path}")

        # ----------------------------------------------------------------------
        # 2. VIRUSES_percentages.tiff (Horizontal bar plot of positivity %)
        # ----------------------------------------------------------------------
        pct_plot_name = f"{fn_prefix}_virus_stats_VIRUSES_percentages.tiff"
        pct_plot_path = os.path.join(out_dir, pct_plot_name)

        v_names_labels = []
        v_positivity_pcts = []
        for v_key in group_v_keys:
            v_info = virus_data[v_key]
            pos_c = len(v_info["samples"])
            ds_k = (v_key[0], v_key[1], v_key[2])
            tot_c = len(total_samples_map[ds_k])
            pct_c = (pos_c / float(tot_c) * 100.0) if tot_c > 0 else 0.0
            lbl = f"{v_info['name']} ({v_key[3]})"
            v_names_labels.append(lbl)
            v_positivity_pcts.append(pct_c)

        if v_names_labels:
            df_pct = pd.DataFrame({'viruses': v_names_labels, 'positivity_pct': v_positivity_pcts})
            df_pct_sorted = df_pct.sort_values(by='positivity_pct', ascending=True)

            num_v = len(df_pct_sorted)
            fig_height = max(8.0, num_v * 0.4 + 2.0)
            fig, ax = plt.subplots(figsize=(10, fig_height), dpi=300)
            bars = ax.barh(df_pct_sorted['viruses'], df_pct_sorted['positivity_pct'], height=0.98, color='#ff98ff')

            ax.set_xlabel('Viral Prevalence (%)', fontsize=16, color='black', labelpad=8)
            ax.set_ylabel('Virus', fontsize=16, color='black', labelpad=8)
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
            ax.set_title(
                f"{dataset} Overall Viral Prevalence (% of Total Cohort Samples)\n"
                f"Phase: {cur_phase_tag} | Project: {cur_proj_tag}",
                fontsize=19, pad=14, color='black'
            )
            max_pct = df_pct_sorted['positivity_pct'].max() if not df_pct_sorted.empty else 1.0
            ax.set_xlim(0, max_pct * 1.22)

            for bar in bars:
                w = bar.get_width()
                ax.annotate(f"{w:.2f}%", xy=(w, bar.get_y() + bar.get_height() / 2),
                            xytext=(6, 0), textcoords="offset points", ha='left', va='center', fontsize=13, color='black')

            ax.grid(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.savefig(pct_plot_path, format="tiff", dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"[SUCCESS] Generated VIRUSES percentages TIFF plot: {pct_plot_path}")

        # ----------------------------------------------------------------------
        # 3. VIRUSES_reads.tiff (Horizontal bar plot of total assigned reads)
        # ----------------------------------------------------------------------
        reads_plot_name = f"{fn_prefix}_virus_stats_VIRUSES_reads.tiff"
        reads_plot_path = os.path.join(out_dir, reads_plot_name)

        v_reads_counts = []
        for v_key in group_v_keys:
            v_info = virus_data[v_key]
            tot_reads_val = sum(v_info["samples"].values())
            v_reads_counts.append(tot_reads_val)

        if v_names_labels:
            df_reads = pd.DataFrame({'viruses': v_names_labels, 'total_reads': v_reads_counts})
            df_reads_sorted = df_reads.sort_values(by='total_reads', ascending=True)

            fig, ax = plt.subplots(figsize=(10, fig_height), dpi=300)
            bars = ax.barh(df_reads_sorted['viruses'], df_reads_sorted['total_reads'], height=0.98, color='#ff98ff')

            ax.set_xlabel('Total Mapped Reads', fontsize=16, color='black', labelpad=8)
            ax.set_ylabel('Virus', fontsize=16, color='black', labelpad=8)
            ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
            ax.set_title(
                f"{dataset} Overall Mapped Read Count Across Cohort\n"
                f"Phase: {cur_phase_tag} | Project: {cur_proj_tag}",
                fontsize=19, pad=14, color='black'
            )
            max_r = df_reads_sorted['total_reads'].max() if not df_reads_sorted.empty else 1.0
            ax.set_xlim(0, max_r * 1.22)

            for bar in bars:
                w = bar.get_width()
                ax.annotate(f"{int(w):,}", xy=(w, bar.get_y() + bar.get_height() / 2),
                            xytext=(6, 0), textcoords="offset points", ha='left', va='center', fontsize=13, color='black')

            ax.grid(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.savefig(reads_plot_path, format="tiff", dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"[SUCCESS] Generated VIRUSES reads TIFF plot: {reads_plot_path}")

        # ----------------------------------------------------------------------
        # 4. Overall Read Counts (2-Panel Figure): {dataset}_{genome_build}_{phase}_{project}_{strategy}_virus_stats_reads.tiff
        # ----------------------------------------------------------------------
        ds_group_key = (cur_phase, cur_proj)
        all_reads_list = dataset_records[ds_group_key]["reads"]

        if all_reads_list:
            all_reads_list = sorted(all_reads_list)
            tot_v_pos = len(all_reads_list)
            reads_overall_name = f"{fn_prefix}_virus_stats_reads.tiff"
            reads_overall_path = os.path.join(out_dir, reads_overall_name)

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5), dpi=300)
            ax1.set_box_aspect(1)
            ax2.set_box_aspect(1)
            overall_color = "#008fbf"

            use_log_a = (args.reads_panel_a_loglog == "on")
            use_log_b = (args.reads_panel_b_log_y == "on")

            u_vals, u_cnts = np.unique(all_reads_list, return_counts=True)
            max_r = max(u_vals) if len(u_vals) > 0 else 1
            min_r = min(u_vals) if len(u_vals) > 0 else 1

            if use_log_a and max_r > args.log_scale_read_cutoff and len(u_vals) > 1:
                bins = np.logspace(np.log10(max(1, min_r)), np.log10(max_r), 30)
                ax1.hist(all_reads_list, bins=bins, color=overall_color, rwidth=0.92, edgecolor='none')
                ax1.set_xscale('log')
                ax1.set_yscale('log')
                ax1.set_xlabel("Mapped Read Count (All Viruses, Log Scale)", fontsize=11, labelpad=8, color='black')
                ax1.set_ylabel("Number of Mapped Occurrences (Log Scale)", fontsize=11, labelpad=8, color='black')
                ax1.set_title("A. Dataset-Wide Read Count Frequency (Log-Log)", fontsize=12, pad=12, color='black')
            else:
                ax1.set_xscale('linear')
                ax1.set_yscale('linear')
                ax1.bar(u_vals, u_cnts, width=0.88, color=overall_color)
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
            bar_w_o = 1.0 if tot_v_pos > 50 else 0.88
            ax2.bar(o_ranks, all_reads_list, width=bar_w_o, color=overall_color, edgecolor='none')
            
            if use_log_b and max_r > args.log_scale_read_cutoff and len(u_vals) > 1:
                ax2.set_yscale('log')
                ax2.set_ylabel("Mapped Read Count (Log Scale)", fontsize=11, labelpad=8, color='black')
            else:
                ax2.set_yscale('linear')
                ax2.set_ylabel("Mapped Read Count", fontsize=11, labelpad=8, color='black')
                ax2.yaxis.set_major_locator(MaxNLocator(integer=True))
                ax2.yaxis.set_major_formatter(FormatStrFormatter('%d'))

            ax2.set_xlabel("Occurrence Rank Across Entire Dataset", fontsize=11, labelpad=8, color='black')
            ax2.set_title("B. Dataset-Wide Ordered Read Counts", fontsize=12, pad=12, color='black')
            ax2.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax2.grid(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)

            fig.suptitle(f"{dataset} Dataset-Wide Read Counts Summary | Phase: {cur_phase_tag} | Project: {cur_proj_tag}\nTotal Viral Detections N = {tot_v_pos}", fontsize=13, y=1.02, color='black')
            fig.text(0.5, 0.005, "Overall dataset-wide read count distribution across all detected viruses", ha='center', fontsize=10.5, style='normal', color='black')
            fig.subplots_adjust(top=0.82, bottom=0.18, left=0.10, right=0.95, wspace=0.30)
            plt.savefig(reads_overall_path, format="tiff", dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"[SUCCESS] Generated overall reads TIFF plot: {reads_overall_path}")

        # ----------------------------------------------------------------------
        # 5. Overall Copy Number (2-Panel Figure): {dataset}_{genome_build}_{phase}_{project}_{strategy}_virus_stats_copy_number.tiff
        # ----------------------------------------------------------------------
        all_cn_list = dataset_records[ds_group_key]["copy_numbers"]

        if all_cn_list:
            all_cn_list = sorted(all_cn_list)
            tot_v_pos = len(all_cn_list)
            cn_overall_name = f"{fn_prefix}_virus_stats_copy_number.tiff"
            cn_overall_path = os.path.join(out_dir, cn_overall_name)

            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5.5), dpi=300)
            ax1.set_box_aspect(1)
            ax2.set_box_aspect(1)
            overall_color = "#008fbf"

            use_log_a = (args.copy_number_panel_a_loglog == "on")
            use_log_b = (args.copy_number_panel_b_log_y == "on")

            u_vals, u_cnts = np.unique(all_cn_list, return_counts=True)
            max_r = max(u_vals) if len(u_vals) > 0 else 1.0
            min_r = min(u_vals) if len(u_vals) > 0 else 0.001

            cn_cutoff = float(getattr(args, "log_scale_copy_number_cutoff", 0.1))

            if use_log_a and max_r > cn_cutoff and len(u_vals) > 1:
                bins = np.logspace(np.log10(max(1e-4, min_r)), np.log10(max_r), 30)
                ax1.hist(all_cn_list, bins=bins, color=overall_color, rwidth=0.92, edgecolor='none')
                ax1.set_xscale('log')
                ax1.set_yscale('log')
                ax1.set_xlabel("Viral Copy Number (Log Scale)", fontsize=11, labelpad=8, color='black')
                ax1.set_ylabel("Number of Mapped Occurrences (Log Scale)", fontsize=11, labelpad=8, color='black')
                ax1.set_title("A. Dataset-Wide Copy Number Frequency (Log-Log)", fontsize=12, pad=12, color='black')
            else:
                ax1.set_xscale('linear')
                ax1.set_yscale('linear')
                ax1.hist(all_cn_list, bins=20, color=overall_color, rwidth=0.92, edgecolor='none')
                ax1.set_xlabel("Viral Copy Number", fontsize=11, labelpad=8, color='black')
                ax1.set_ylabel("Number of Mapped Occurrences", fontsize=11, labelpad=8, color='black')
                ax1.set_title("A. Dataset-Wide Copy Number Frequency", fontsize=12, pad=12, color='black')

            ax1.grid(False)
            ax1.spines['top'].set_visible(False)
            ax1.spines['right'].set_visible(False)

            o_ranks = np.arange(1, tot_v_pos + 1)
            bar_w_o = 1.0 if tot_v_pos > 50 else 0.88
            ax2.bar(o_ranks, all_cn_list, width=bar_w_o, color=overall_color, edgecolor='none')
            
            if use_log_b and max_r > cn_cutoff and len(u_vals) > 1:
                ax2.set_yscale('log')
                ax2.set_ylabel("Viral Copy Number (Log Scale)", fontsize=11, labelpad=8, color='black')
            else:
                ax2.set_yscale('linear')
                ax2.set_ylabel("Viral Copy Number", fontsize=11, labelpad=8, color='black')

            ax2.set_xlabel("Occurrence Rank Across Entire Dataset", fontsize=11, labelpad=8, color='black')
            ax2.set_title("B. Dataset-Wide Ordered Copy Numbers", fontsize=12, pad=12, color='black')
            ax2.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
            ax2.xaxis.set_major_formatter(FormatStrFormatter('%d'))
            ax2.grid(False)
            ax2.spines['top'].set_visible(False)
            ax2.spines['right'].set_visible(False)

            fig.suptitle(f"{dataset} Dataset-Wide Copy Number Summary | Phase: {cur_phase_tag} | Project: {cur_proj_tag}\nTotal Viral Detections N = {tot_v_pos}", fontsize=13, y=1.02, color='black')
            fig.text(0.5, 0.005, "Overall dataset-wide copy number distribution across all detected viruses", ha='center', fontsize=10.5, style='normal', color='black')
            fig.subplots_adjust(top=0.82, bottom=0.18, left=0.10, right=0.95, wspace=0.30)
            plt.savefig(cn_overall_path, format="tiff", dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"[SUCCESS] Generated overall copy number TIFF plot: {cn_overall_path}")

    # --------------------------------------------------------------------------
    # 6. Heatmaps Module Integration (Outputs TIFF heatmaps in stats out_dir)
    # --------------------------------------------------------------------------
    rc_switch = str(getattr(args, "heatmap_read_counts", "off")).lower()
    cn_switch = str(getattr(args, "heatmap_copy_number", "off")).lower()
    hm_strategy = getattr(args, "target_heatmap_strategy", "clean_flags")

    if rc_switch == "on" or cn_switch == "on":
        import subprocess
        heatmap_script_path = os.path.join(os.path.dirname(__file__), "generate_heatmap.py")
        if not os.path.exists(heatmap_script_path):
            heatmap_script_path = "scripts/generate_heatmap.py"

        if rc_switch == "on":
            cmd_rc = (
                f"python3 {heatmap_script_path} --input-report \"{input_report_path}\" "
                f"--out-dir \"{out_dir}\" --dataset \"{dataset}\" --genome-build \"{genome_build}\" "
                f"--phase \"{args.target_phase or ''}\" --project \"{args.target_project or ''}\" "
                f"--strategy \"{hm_strategy}\" --value-type \"read_counts\""
            )
            subprocess.run(cmd_rc, shell=True, check=True)

        if cn_switch == "on":
            cmd_cn = (
                f"python3 {heatmap_script_path} --input-report \"{input_report_path}\" "
                f"--out-dir \"{out_dir}\" --dataset \"{dataset}\" --genome-build \"{genome_build}\" "
                f"--phase \"{args.target_phase or ''}\" --project \"{args.target_project or ''}\" "
                f"--strategy \"{hm_strategy}\" --value-type \"copy_number\""
            )
            subprocess.run(cmd_cn, shell=True, check=True)


def main():
    args = parse_args()
    generate_stats(args)
    from datetime import datetime
    end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n[COMPLETED] Execution finished at: {end_time}")


if __name__ == "__main__":
    main()

