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
    parser.add_argument("--heatmap-min-reads-cutoff", type=int, default=3,
                        help="Minimum mapped reads threshold for heatmap inclusion (default: 3)")
    parser.add_argument("--heatmap-stratify-by", default="",
                        help="Demographics column to split heatmaps by (e.g., sex, phenotype)")
    parser.add_argument("--group-level", default="species", choices=["none", "species", "genus", "family", "realm"],
                        help="Taxonomic grouping level for aggregated summary TSV and bar plots (default: species)")
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
        'sample_id': [c for c in df_raw.columns if 'sample' in c][0] if any('sample' in c for c in df_raw.columns) else df_raw.columns[0],
        'virus_accession': [c for c in df_raw.columns if 'accession' in c][0] if any('accession' in c for c in df_raw.columns) else df_raw.columns[1],
        'virus_mapped_reads': [c for c in df_raw.columns if 'mapped_reads' in c][0] if any('mapped_reads' in c for c in df_raw.columns) else df_raw.columns[3],
        'viral_copy_number': [c for c in df_raw.columns if 'copy_number' in c][0] if any('copy_number' in c for c in df_raw.columns) else df_raw.columns[8],
        'virus_name_sanitized': [c for c in df_raw.columns if 'name_sanitized' in c or 'virus_name' in c][0] if any('name' in c for c in df_raw.columns) else df_raw.columns[9],
        'source_file': [c for c in df_raw.columns if 'source_file' in c or 'source' in c][0] if any('source' in c for c in df_raw.columns) else df_raw.columns[10],
        'phase': [c for c in df_raw.columns if 'phase' in c][0] if any('phase' in c for c in df_raw.columns) else df_raw.columns[11],
        'project': [c for c in df_raw.columns if 'project' in c][0] if any('project' in c for c in df_raw.columns) else df_raw.columns[12],
        'species_taxid': [c for c in df_raw.columns if 'species_taxid' in c][0] if any('species_taxid' in c for c in df_raw.columns) else '',
        'species_name': [c for c in df_raw.columns if 'species_name' in c][0] if any('species_name' in c for c in df_raw.columns) else '',
        'genus_name': [c for c in df_raw.columns if 'genus_name' in c][0] if any('genus_name' in c for c in df_raw.columns) else '',
        'family_name': [c for c in df_raw.columns if 'family_name' in c][0] if any('family_name' in c for c in df_raw.columns) else '',
        'realm_name': [c for c in df_raw.columns if 'realm_name' in c][0] if any('realm_name' in c for c in df_raw.columns) else ''
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

    df['species_taxid'] = df_raw[col_map['species_taxid']].fillna('Unknown').astype(str).str.strip() if col_map['species_taxid'] else 'Unknown'
    df['species_name'] = df_raw[col_map['species_name']].fillna('Unknown').astype(str).str.strip() if col_map['species_name'] else 'Unknown'
    df['genus_name'] = df_raw[col_map['genus_name']].fillna('Unknown').astype(str).str.strip() if col_map['genus_name'] else 'Unknown'
    df['family_name'] = df_raw[col_map['family_name']].fillna('Unknown').astype(str).str.strip() if col_map['family_name'] else 'Unknown'
    df['realm_name'] = df_raw[col_map['realm_name']].fillna('Unknown').astype(str).str.strip() if col_map['realm_name'] else 'Unknown'

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
        "species_taxid": "Unknown",
        "species_name": "Unknown",
        "genus_name": "Unknown",
        "family_name": "Unknown",
        "realm_name": "Unknown",
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
                virus_data[v_key]["species_taxid"] = row.species_taxid
                virus_data[v_key]["species_name"] = row.species_name
                virus_data[v_key]["genus_name"] = row.genus_name
                virus_data[v_key]["family_name"] = row.family_name
                virus_data[v_key]["realm_name"] = row.realm_name
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

        strat_tag = get_short_strategy(group_v_keys[0][2]) if group_v_keys else "clean_flags"

        p_clean = str(cur_phase)[5:] if str(cur_phase).lower().startswith("phase") else str(cur_phase)
        p_tag = "ALL" if is_phase_empty or not p_clean else p_clean
        prj_tag = "ALL" if is_proj_empty else str(cur_proj)
        cur_phase_tag = p_tag
        cur_proj_tag = prj_tag

        fn_prefix = f"{dataset}_{genome_build}_{p_tag}_{prj_tag}_{strat_tag}"

        # ----------------------------------------------------------------------
        # 1. Generate VIRUSES TSV & VIRUSES_classes.tsv
        # ----------------------------------------------------------------------
        viruses_tsv_name = f"{fn_prefix}_virus_stats_VIRUSES.tsv"
        viruses_tsv_path = os.path.join(out_dir, viruses_tsv_name)
        classes_tsv_name = f"{fn_prefix}_virus_stats_VIRUSES_classes.tsv"
        classes_tsv_path = os.path.join(out_dir, classes_tsv_name)

        viruses_header = [
            "Phase", "Project", "Strategy", "Virus_Accession", "Virus_Name_Sanitized",
            "Species_TaxID", "Species_Name", "Genus_Name", "Family_Name", "Realm_Name",
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

        # Pre-sort group_v_keys hierarchically by taxonomy: family -> genus -> species_taxid -> virus_accession
        sorted_group_v_keys = sorted(
            group_v_keys,
            key=lambda x: (
                virus_data[x]["family_name"],
                virus_data[x]["genus_name"],
                virus_data[x]["species_taxid"],
                x[3]
            )
        )

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
                v_info["species_taxid"], v_info["species_name"], v_info["genus_name"],
                v_info["family_name"], v_info["realm_name"],
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
        # 1b. Grouped Summary TSV & Group Counts Bar Plot (--group-level)
        # ----------------------------------------------------------------------
        grp_lvl = args.group_level.strip().lower()
        if grp_lvl != "none":
            grp_tsv_name = f"{fn_prefix}_virus_stats_grouped_{grp_lvl}.tsv"
            grp_tsv_path = os.path.join(out_dir, grp_tsv_name)
            
            grp_counts = collections.defaultdict(lambda: {"pos_samples": set(), "total_reads": 0})
            for v_key in sorted_group_v_keys:
                v_info = virus_data[v_key]
                if grp_lvl == "species":
                    g_key = f"{v_info['species_taxid']}|{v_info['species_name']}"
                elif grp_lvl == "genus":
                    g_key = v_info["genus_name"]
                elif grp_lvl == "family":
                    g_key = v_info["family_name"]
                elif grp_lvl == "realm":
                    g_key = v_info["realm_name"]
                else:
                    g_key = f"{v_info['species_taxid']}|{v_info['species_name']}"

                grp_counts[g_key]["pos_samples"].update(v_info["samples"].keys())
                grp_counts[g_key]["total_reads"] += sum(v_info["samples"].values())

            grp_rows = []
            grp_header = [grp_lvl.capitalize() + "_Group", "Positive_Samples_Count", "Total_Reads_Assigned"]
            for g_k, g_v in sorted(grp_counts.items(), key=lambda x: len(x[1]["pos_samples"]), reverse=True):
                grp_rows.append([g_k, str(len(g_v["pos_samples"])), str(g_v["total_reads"])])

            with open(grp_tsv_path, "w") as f:
                f.write("\t".join(grp_header) + "\n")
                for r in grp_rows:
                    f.write("\t".join(r) + "\n")
            print(f"[SUCCESS] Written grouped summary TSV: {grp_tsv_path}")

            # Render Group Counts Bar Plot: {dataset}_{genome_build}_{phase}_{project}_stats_{group_level}_counts.tiff
            grp_plot_name = f"{fn_prefix}_stats_{grp_lvl}_counts.tiff"
            grp_plot_path = os.path.join(out_dir, grp_plot_name)

            top_grp_items = sorted(grp_counts.items(), key=lambda x: len(x[1]["pos_samples"]), reverse=False)
            g_labels = [x[0] for x in top_grp_items]
            g_pos_counts = [len(x[1]["pos_samples"]) for x in top_grp_items]

            if g_labels:
                fig_h = max(7.0, len(g_labels) * 0.4 + 2.0)
                fig, ax = plt.subplots(figsize=(10, fig_h), dpi=300)
                bars = ax.barh(g_labels, g_pos_counts, height=0.88, color="#008fbf")
                ax.set_xlabel("Positive Sample Detections", fontsize=14, labelpad=8)
                ax.set_ylabel(f"Taxonomic Group ({grp_lvl.capitalize()})", fontsize=14, labelpad=8)
                ax.tick_params(axis='x', labelsize=12)
                ax.tick_params(axis='y', labelsize=12)
                ax.set_title(f"{dataset} Cohort Positivity by Taxonomic {grp_lvl.capitalize()}\nPhase: {cur_phase_tag} | Project: {cur_proj_tag}", fontsize=16, pad=12)
                
                max_g = max(g_pos_counts) if g_pos_counts else 1
                ax.set_xlim(0, max_g * 1.18)
                for bar in bars:
                    w = bar.get_width()
                    ax.annotate(f"{int(w)}", xy=(w, bar.get_y() + bar.get_height() / 2), xytext=(5, 0), textcoords="offset points", ha='left', va='center', fontsize=12)
                
                ax.grid(False)
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                plt.savefig(grp_plot_path, format="tiff", dpi=300, bbox_inches="tight")
                plt.close(fig)
                print(f"[SUCCESS] Generated grouped counts bar plot: {grp_plot_path}")

        # ----------------------------------------------------------------------
        # 1c. VIRUSES_classes.tiff (Vertical bar plot of 4 classification categories ordered highest to lowest)
        # ----------------------------------------------------------------------
        classes_plot_name = f"{fn_prefix}_virus_stats_VIRUSES_classes.tiff"
        classes_plot_path = os.path.join(out_dir, classes_plot_name)

        sorted_classes = class_counts.most_common()
        cls_names = [c[0] for c in sorted_classes]
        cls_counts = [c[1] for c in sorted_classes]

        fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
        bars = ax.bar(cls_names, cls_counts, width=0.95, color="#008fbf", edgecolor='none')
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
        # Format: Species_TaxID|Species_Name|Virus_Accession|Virus_Name_Sanitized
        # Ordered by taxonomy: family -> genus -> species_taxid -> virus_accession
        # ----------------------------------------------------------------------
        pct_plot_name = f"{fn_prefix}_virus_stats_VIRUSES_percentages.tiff"
        pct_plot_path = os.path.join(out_dir, pct_plot_name)

        v_names_labels = []
        v_positivity_pcts = []
        v_reads_counts = []
        for v_key in sorted_group_v_keys:
            v_info = virus_data[v_key]
            pos_c = len(v_info["samples"])
            ds_k = (v_key[0], v_key[1], v_key[2])
            tot_c = len(total_samples_map[ds_k])
            pct_c = (pos_c / float(tot_c) * 100.0) if tot_c > 0 else 0.0
            lbl = f"{v_info['species_taxid']}|{v_info['species_name']}|{v_key[3]}|{v_info['name']}"
            v_names_labels.append(lbl)
            v_positivity_pcts.append(pct_c)
            v_reads_counts.append(sum(v_info["samples"].values()))

        if v_names_labels:
            df_pct = pd.DataFrame({'viruses': v_names_labels, 'positivity_pct': v_positivity_pcts})

            num_v = len(df_pct)
            fig_height = max(8.0, num_v * 0.45 + 2.0)
            fig, ax = plt.subplots(figsize=(12, fig_height), dpi=300)
            bars = ax.barh(df_pct['viruses'], df_pct['positivity_pct'], height=0.88, color='#ff98ff')

            ax.set_xlabel('Viral Prevalence (%)', fontsize=16, color='black', labelpad=8)
            ax.set_ylabel('Viral Reference (TaxID|Species|Accession|Name)', fontsize=16, color='black', labelpad=8)
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=11)
            ax.set_title(
                f"{dataset} Overall Viral Prevalence (% of Total Cohort Samples)\n"
                f"Phase: {cur_phase_tag} | Project: {cur_proj_tag}",
                fontsize=18, pad=14, color='black'
            )
            max_pct = df_pct['positivity_pct'].max() if not df_pct.empty else 1.0
            ax.set_xlim(0, max_pct * 1.22)

            for bar in bars:
                w = bar.get_width()
                ax.annotate(f"{w:.2f}%", xy=(w, bar.get_y() + bar.get_height() / 2),
                            xytext=(6, 0), textcoords="offset points", ha='left', va='center', fontsize=11, color='black')

            ax.grid(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.savefig(pct_plot_path, format="tiff", dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"[SUCCESS] Generated VIRUSES percentages TIFF plot: {pct_plot_path}")

        # ----------------------------------------------------------------------
        # 3. VIRUSES_reads.tiff (Horizontal bar plot of total assigned reads)
        # Format: Species_TaxID|Species_Name|Virus_Accession|Virus_Name_Sanitized
        # Ordered by taxonomy: family -> genus -> species_taxid -> virus_accession
        # ----------------------------------------------------------------------
        reads_plot_name = f"{fn_prefix}_virus_stats_VIRUSES_reads.tiff"
        reads_plot_path = os.path.join(out_dir, reads_plot_name)

        if v_names_labels:
            df_reads = pd.DataFrame({'viruses': v_names_labels, 'total_reads': v_reads_counts})

            fig, ax = plt.subplots(figsize=(12, fig_height), dpi=300)
            bars = ax.barh(df_reads['viruses'], df_reads['total_reads'], height=0.88, color='#ff98ff')

            ax.set_xlabel('Total Mapped Reads', fontsize=16, color='black', labelpad=8)
            ax.set_ylabel('Viral Reference (TaxID|Species|Accession|Name)', fontsize=16, color='black', labelpad=8)
            ax.xaxis.set_major_formatter(matplotlib.ticker.StrMethodFormatter('{x:,.0f}'))
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=11)
            ax.set_title(
                f"{dataset} Overall Mapped Read Count Across Cohort\n"
                f"Phase: {cur_phase_tag} | Project: {cur_proj_tag}",
                fontsize=18, pad=14, color='black'
            )
            max_r = df_reads['total_reads'].max() if not df_reads.empty else 1.0
            ax.set_xlim(0, max_r * 1.22)

            for bar in bars:
                w = bar.get_width()
                ax.annotate(f"{int(w):,}", xy=(w, bar.get_y() + bar.get_height() / 2),
                            xytext=(6, 0), textcoords="offset points", ha='left', va='center', fontsize=11, color='black')

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
    # 5. Taxonomy & Demographics Plots Integration (from OLD/taxonomy_demographics_plots.py)
    # --------------------------------------------------------------------------
    import matplotlib.patches as mpatches
    import matplotlib.ticker as mtick

    PARENT_GROUP = {
        "Accession": "species_name",
        "Species":   "genus_name",
        "Genus":     "family_name",
        "Family":    "order_name",
        "Order":     "realm_name",
        "Realm":     "realm_name",
    }

    PALETTE = [
        "#0072B2", "#E69F00", "#009E73", "#D55E00", "#7570B3", "#56B4E9", "#E7298A", "#66A61E",
        "#1B9E77", "#CC79A7", "#A6761D", "#4B0082", "#FC8D62", "#8DA0CB", "#666666", "#E78AC3",
        "#A6D854", "#984EA3", "#377EB8", "#FFD92F", "#E5C494", "#B3DE69", "#BC80BD", "#FB8072",
        "#80B1D3", "#B3B3B3"
    ]

    # Run taxonomy level and demographics plots for filtered df
    df_plot_filt = df[df["virus_accession"].notna() & (df["virus_accession"] != "None") & (df["virus_accession"] != "")].copy()
    
    # Rescue unclassified taxonomy names
    mask_pap = df_plot_filt["virus_name_sanitized"].str.contains("papilloma", case=False, na=False)
    df_plot_filt.loc[mask_pap & (df_plot_filt["family_name"] == "Unknown"), "family_name"] = "Papillomaviridae"
    df_plot_filt.loc[mask_pap & (df_plot_filt["genus_name"] == "Unknown"), "genus_name"] = "Betapapillomavirus"
    df_plot_filt.loc[mask_pap & (df_plot_filt["order_name"] == "Unknown"), "order_name"] = "Zurhausenvirales"
    df_plot_filt.loc[mask_pap & (df_plot_filt["realm_name"] == "Unknown"), "realm_name"] = "Floreoviria"

    mask_herpes = df_plot_filt["virus_name_sanitized"].str.contains("herpes", case=False, na=False)
    df_plot_filt.loc[mask_herpes & (df_plot_filt["family_name"] == "Unknown"), "family_name"] = "Orthoherpesviridae"
    df_plot_filt.loc[mask_herpes & (df_plot_filt["order_name"] == "Unknown"), "order_name"] = "Herpesvirales"
    df_plot_filt.loc[mask_herpes & (df_plot_filt["realm_name"] == "Unknown"), "realm_name"] = "Duplodnaviria"

    mask_ttv = df_plot_filt["virus_name_sanitized"].str.contains("torque teno|anellovir", case=False, na=False)
    df_plot_filt.loc[mask_ttv & (df_plot_filt["family_name"] == "Unknown"), "family_name"] = "Anelloviridae"

    if not df_plot_filt.empty:
        n_tot = df_plot_filt["sample_id"].nunique()

        # 5a. Taxonomy Horizontal Bar Plots at 6 levels
        tax_levels = [
            ("virus_accession",  "virus_name_sanitized", "Accession",  "accession"),
            ("species_taxid",    "species_name",         "Species",    "species"),
            ("genus_taxid",      "genus_name",           "Genus",      "genus"),
            ("family_taxid",     "family_name",          "Family",     "family"),
            ("order_taxid",      "order_name",           "Order",      "order"),
            ("realm_taxid",      "realm_name",           "Realm",      "realm"),
        ]

        for taxid_c, name_c, lvl_lbl, lvl_fname in tax_levels:
            if taxid_c not in df_plot_filt.columns or name_c not in df_plot_filt.columns:
                continue

            p_data = df_plot_filt.copy()
            if lvl_lbl != "Accession":
                p_data = p_data[p_data[name_c] != "Unknown"]

            parent_c = PARENT_GROUP[lvl_lbl]
            above_ranks = []
            if lvl_lbl == "Accession":
                above_ranks = ["family_name", "genus_name", "species_name"]
            elif lvl_lbl == "Species":
                above_ranks = ["family_name", "genus_name"]
            elif lvl_lbl == "Genus":
                above_ranks = ["family_name"]

            grp_c = [taxid_c, name_c]
            for c in above_ranks:
                if c not in grp_c and c in p_data.columns:
                    grp_c.append(c)
            if parent_c not in grp_c and parent_c in p_data.columns:
                grp_c.append(parent_c)

            if p_data.empty:
                continue

            grp = p_data.groupby(grp_c).agg(
                n_positive_samples=("sample_id", "nunique"),
                total_reads=("virus_mapped_reads", "sum")
            ).reset_index()

            grp["pct_positive"] = grp["n_positive_samples"] / n_tot * 100.0
            grp["label"] = grp[taxid_c].astype(str) + " – " + grp[name_c].astype(str)

            s_cols = []
            s_asc = []
            for c in above_ranks:
                if c in p_data.columns:
                    c_tot = p_data.groupby(c)["sample_id"].nunique().to_dict()
                    grp[c + "_tot"] = grp[c].map(c_tot)
                    grp[c + "_unk"] = grp[c] == "Unknown"
                    s_cols.extend([c + "_unk", c + "_tot", c])
                    s_asc.extend([True, False, True])

            s_cols.extend(["n_positive_samples", "total_reads"])
            s_asc.extend([False, False])
            grp = grp.sort_values(by=s_cols, ascending=s_asc).reset_index(drop=True)

            uniq_grps = list(dict.fromkeys(grp[parent_c])) if parent_c in grp.columns else ["Unknown"]
            color_m = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(uniq_grps)}
            bar_cols = [color_m.get(g, PALETTE[0]) for g in (grp[parent_c] if parent_c in grp.columns else ["Unknown"]*len(grp))]

            fig_h = max(3.5, len(grp) * 0.35)
            fig, ax = plt.subplots(figsize=(10.5, fig_h))
            bars = ax.barh(grp["label"], grp["pct_positive"], color=bar_cols, edgecolor="none", height=0.90)
            ax.invert_yaxis()

            max_v = grp["pct_positive"].max() if len(grp) else 1
            for bar, reads in zip(bars, grp["total_reads"]):
                w = bar.get_width()
                if w > max_v * 0.25:
                    ax.text(w - max_v * 0.015, bar.get_y() + bar.get_height()/2, f"{int(reads):,} reads", va="center", ha="right", fontsize=8.5, color="white")
                else:
                    ax.text(w + max_v * 0.015, bar.get_y() + bar.get_height()/2, f"{int(reads):,} reads", va="center", ha="left", fontsize=8.5, color="black")

            handles = [mpatches.Patch(facecolor=color_m[g], edgecolor="none", label=g) for g in uniq_grps]
            ax.legend(handles=handles, title=parent_c.replace("_name", "").title(), fontsize=8, title_fontsize=8.5, bbox_to_anchor=(1.01, 0.5), loc="center left", frameon=False)
            ax.set_xlabel(f"% positive samples (n = {n_tot})", fontsize=10)
            ax.set_title(f"Distribution of Viruses – {lvl_lbl} Level\nPhase: {cur_phase_tag} | Project: {cur_proj_tag}", fontsize=10.5)
            ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
            ax.tick_params(axis="x", labelsize=9)
            ax.tick_params(axis="y", labelsize=8.5)
            fig.text(0.5, 0.005, f"source: {strat_tag}", ha="center", fontsize=9.5, color="black")
            plt.tight_layout(rect=[0, 0.04, 1, 1])

            t_out_name = f"{fn_prefix}_taxonomy_{lvl_fname}.tiff"
            t_out_path = os.path.join(out_dir, t_out_name)
            fig.savefig(t_out_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
            plt.close(fig)
            print(f"[SUCCESS] Saved taxonomy bar plot: {t_out_path}")

        # 5b. Accession Plot with Family Overlay
        above = ["family_name", "genus_name", "species_name"]
        grp_c = ["virus_accession", "virus_name_sanitized"] + [c for c in above if c in df_plot_filt.columns]
        grp = df_plot_filt.groupby(grp_c).agg(n_positive_samples=("sample_id", "nunique"), total_reads=("virus_mapped_reads", "sum")).reset_index()
        grp["pct_positive"] = grp["n_positive_samples"] / n_tot * 100.0
        grp["label"] = grp["virus_accession"].astype(str) + " – " + grp["virus_name_sanitized"].astype(str)

        s_cols, s_asc = [], []
        for c in above:
            if c in df_plot_filt.columns:
                c_tot = df_plot_filt.groupby(c)["sample_id"].nunique().to_dict()
                grp[c + "_tot"] = grp[c].map(c_tot)
                grp[c + "_unk"] = grp[c] == "Unknown"
                s_cols.extend([c + "_unk", c + "_tot", c])
                s_asc.extend([True, False, True])
        s_cols.extend(["n_positive_samples", "total_reads"])
        s_asc.extend([False, False])
        grp = grp.sort_values(by=s_cols, ascending=s_asc).reset_index(drop=True)

        if not grp.empty and "genus_name" in grp.columns and "family_name" in grp.columns:
            genus_order = list(dict.fromkeys(grp["genus_name"]))
            genus_color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(genus_order)}
            acc_bar_colors = [genus_color_map[g] for g in grp["genus_name"]]

            fam_order = list(dict.fromkeys(grp["family_name"]))
            fam_samples = df_plot_filt.groupby("family_name")["sample_id"].nunique().to_dict()
            fam_reads = df_plot_filt.groupby("family_name")["virus_mapped_reads"].sum().to_dict()

            fam_palette = ["#C6D8EF", "#FEE8C8", "#E7D4E8", "#FDDBC7", "#E0ECF4", "#E5F5E0"]
            fam_tint_map = {f: fam_palette[i % len(fam_palette)] for i, f in enumerate(fam_order)}

            fig_h = max(4.8, len(grp) * 0.36)
            fig, ax = plt.subplots(figsize=(11.5, fig_h))
            y_positions = list(range(len(grp)))

            for fam in fam_order:
                idx_list = grp.index[grp["family_name"] == fam].tolist()
                if not idx_list:
                    continue
                y_min = min(idx_list) - 0.46
                y_max = max(idx_list) + 0.46
                pct = fam_samples[fam] / n_tot * 100.0
                tint = fam_tint_map[fam]
                rect = plt.Rectangle((0, y_min), pct, y_max - y_min, facecolor=tint, alpha=0.6, edgecolor="none", zorder=1)
                ax.add_patch(rect)
                fam_label = f"Family {fam}: {fam_samples[fam]} samples ({pct:.1f}%), {int(fam_reads[fam]):,} reads"
                x_pos = pct + 1.0 if pct > 10 else 13.0
                ax.text(x_pos, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=8.5, color="#1A252F", zorder=2)

            bars = ax.barh(y_positions, grp["pct_positive"], height=0.88, color=acc_bar_colors, edgecolor="none", zorder=3)
            ax.set_yticks(y_positions)
            ax.set_yticklabels(grp["label"], fontsize=8.5)
            ax.invert_yaxis()

            for bar, reads in zip(bars, grp["total_reads"]):
                w = bar.get_width()
                if w > 12:
                    ax.text(w - 1, bar.get_y() + bar.get_height() / 2, f"{int(reads):,} reads", va="center", ha="right", fontsize=8.5, color="white", zorder=4)
                else:
                    ax.text(w + 0.25, bar.get_y() + bar.get_height() / 2, f"{int(reads):,} reads", va="center", ha="left", fontsize=8.5, color="black", zorder=4)

            handles = [mpatches.Patch(facecolor=genus_color_map[g], edgecolor="none", label=g) for g in genus_order]
            ax.legend(handles=handles, title="Genus", fontsize=8, title_fontsize=8.5, bbox_to_anchor=(1.01, 0.5), loc="center left", frameon=False)
            ax.set_xlabel(f"% positive samples (n = {n_tot})", fontsize=10)
            ax.set_title(f"Distribution of Viruses – Accession Level with Family Overlay\nPhase: {cur_phase_tag} | Project: {cur_proj_tag}", fontsize=10.5)
            ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
            ax.tick_params(axis="x", labelsize=9)
            ax.tick_params(axis="y", labelsize=8.5)
            max_v_fam = max(fam_samples.values()) / n_tot * 100.0 if fam_samples else 10.0
            ax.set_xlim(0, max_v_fam * 1.55)
            fig.text(0.5, 0.005, f"source: {strat_tag}", ha="center", fontsize=9.5, color="black")
            plt.tight_layout(rect=[0, 0.04, 1, 1])

            overlay_out_name = f"{fn_prefix}_taxonomy_accession_family_overlay.tiff"
            overlay_out_path = os.path.join(out_dir, overlay_out_name)
            fig.savefig(overlay_out_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
            plt.close(fig)
            print(f"[SUCCESS] Saved accession family overlay plot: {overlay_out_path}")

        # 5c. Demographics Plots (Histograms & Vertical Bar Charts)
        demo_cols = [
            ("sex",               "Sex",        "sex",       False),
            ("age",               "Age",        "age",       True),
            ("tissue",            "Tissue",     "tissue",    False),
            ("cell",              "Cell Type",  "cell",      False),
            ("race",              "Race",       "race",      False),
            ("phenotype",         "Phenotype",  "phenotype", False),
            ("sample_read_depth", "Coverage",   "coverage",  True),
        ]

        sample_df = df_plot_filt.drop_duplicates(subset="sample_id")

        for col, label, var_fname, is_num in demo_cols:
            matching_cols = [c for c in sample_df.columns if c.lower().strip() == col]
            if not matching_cols:
                continue

            col_name = matching_cols[0]
            fig, ax = plt.subplots(figsize=(5, 5))

            if is_num:
                vals = pd.to_numeric(sample_df[col_name], errors="coerce").dropna()
                if vals.empty or vals.nunique() <= 1:
                    plt.close(fig)
                    continue
                n_bins = min(20, max(8, int(np.sqrt(len(vals)))))
                ax.hist(vals, bins=n_bins, color="#4472C4", edgecolor="none", rwidth=0.85)
                ax.set_xlabel(label, fontsize=10)
                ax.set_ylabel("Number of samples", fontsize=10)
                ax.tick_params(axis="x", labelsize=9)
            else:
                cat_counts = sample_df[col_name].fillna("Unknown").value_counts()
                if len(cat_counts) <= 1:
                    plt.close(fig)
                    continue
                cat_counts = cat_counts.sort_values(ascending=False)
                x_labels = [str(x).replace("_", " ") for x in cat_counts.index]
                n_cats = len(cat_counts)

                if n_cats == 2:
                    x_pos = [-0.35, 0.35]
                    bars = ax.bar(x_pos, cat_counts.values, color="#4472C4", edgecolor="none", width=0.52)
                    ax.set_xticks(x_pos)
                    ax.set_xticklabels(x_labels, fontsize=9)
                    ax.set_xlim(-1.0, 1.0)
                else:
                    bars = ax.bar(range(n_cats), cat_counts.values, color="#4472C4", edgecolor="none", width=0.75)
                    ax.set_xticks(range(n_cats))
                    ax.set_xticklabels(x_labels, fontsize=9)

                max_v_cat = cat_counts.max() if len(cat_counts) else 1
                for bar, cnt in zip(bars, cat_counts.values):
                    h = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width() / 2, h + max_v_cat * 0.015, f"{cnt}", ha="center", va="bottom", fontsize=8.5, color="black")

                ax.set_ylabel("Number of samples", fontsize=10)
                ax.set_ylim(0, max_v_cat * 1.15)
                ax.tick_params(axis="x", labelsize=9)
                if any(len(lbl) > 8 for lbl in x_labels):
                    plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")

            ax.tick_params(axis="y", labelsize=9)
            ax.set_title(f"Demographics – {label}\nPhase: {cur_phase_tag} | Project: {cur_proj_tag}", fontsize=10.5)
            fig.text(0.5, 0.015, f"source: {strat_tag}", ha="center", fontsize=9.5, color="black")
            plt.tight_layout(rect=[0, 0.05, 1, 1])

            demo_out_name = f"{fn_prefix}_demo_{var_fname}.tiff"
            demo_out_path = os.path.join(out_dir, demo_out_name)
            fig.savefig(demo_out_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
            plt.close(fig)
            print(f"[SUCCESS] Saved demographics plot: {demo_out_path}")

    # --------------------------------------------------------------------------
    # 6. Heatmaps Module Integration (Outputs TIFF heatmaps in stats out_dir)
    # --------------------------------------------------------------------------
    rc_switch = str(getattr(args, "heatmap_read_counts", "off")).lower()
    cn_switch = str(getattr(args, "heatmap_copy_number", "off")).lower()
    hm_min_reads = getattr(args, "heatmap_min_reads_cutoff", 3)
    hm_stratify = getattr(args, "heatmap_stratify_by", "")

    # Derive heatmap strategy from primary target strategy
    primary_strat = args.strategies[0] if args.strategies else "clean_flags"
    hm_strategy = get_short_strategy(primary_strat)

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
                f"--strategy \"{hm_strategy}\" --value-type \"read_counts\" --min-reads-cutoff {hm_min_reads} "
                f"--stratify-by \"{hm_stratify}\""
            )
            subprocess.run(cmd_rc, shell=True, check=True)

        if cn_switch == "on":
            cmd_cn = (
                f"python3 {heatmap_script_path} --input-report \"{input_report_path}\" "
                f"--out-dir \"{out_dir}\" --dataset \"{dataset}\" --genome-build \"{genome_build}\" "
                f"--phase \"{args.target_phase or ''}\" --project \"{args.target_project or ''}\" "
                f"--strategy \"{hm_strategy}\" --value-type \"copy_number\" --min-reads-cutoff {hm_min_reads} "
                f"--stratify-by \"{hm_stratify}\""
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

