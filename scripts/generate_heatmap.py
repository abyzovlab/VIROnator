#!/usr/bin/env python3
"""
generate_heatmap.py
Generates hierarchical clustering heatmaps of viral read counts and/or viral copy numbers
from consolidated master viral report TSV.

Preserves exact formatting & plotting parameters from OLD/06_heatmap.py:
- np.log10(df + 1e-7)
- average linkage clustering on log_df.T
- seaborn heatmap with cmap='viridis', linewidths=0.5
- Title / Label font sizes: Title (70), Labels (120), xticks (100, rot 90), yticks (100, rot -360), colorbar (120, rot 90)
"""

import argparse
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list


def parse_args():
    parser = argparse.ArgumentParser(description="Generate read count and copy number heatmaps from master viral report.")
    parser.add_argument("--input-report", required=True, help="Path to master viral report TSV file")
    parser.add_argument("--out-dir", required=True, help="Output directory for heatmap PNG files")
    parser.add_argument("--dataset", default="DATASET", help="Dataset identifier")
    parser.add_argument("--genome-build", default="hg38", help="Genome build identifier")
    parser.add_argument("--phase", default="", help="Target phase filter from config")
    parser.add_argument("--project", default="", help="Target project filter from config")
    parser.add_argument("--strategy", default="clean_flags", help="Target CRAM strategy string/keyword")
    parser.add_argument("--value-type", choices=["read_counts", "copy_number"], required=True,
                        help="Value type to plot: 'read_counts' or 'copy_number'")
    parser.add_argument("--min-reads-cutoff", type=int, default=3,
                        help="Minimum mapped reads threshold for sample inclusion in heatmap (default: 3)")
    parser.add_argument("--group-level", default="none", choices=["none", "species", "genus", "family", "realm"],
                        help="Taxonomic grouping level for aggregated heatmap (default: none)")
    return parser.parse_args()


def get_short_strategy(fname):
    """Maps long CRAM filename to concise strategy label or returns given string if already concise."""
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
    return fname.strip()


def render_and_save_heatmap(pivot_df, out_path, title_text, value_type, phase, project, strategy):
    """Helper function to log-transform, hierarchically cluster, and save a heatmap."""
    if pivot_df.empty or pivot_df.shape[1] == 0:
        return

    # Log10 transform positive values; zero/negative -> NaN
    log_df = pivot_df.replace(0.0, np.nan)
    log_df = np.log10(log_df)

    num_samples = log_df.shape[0]
    num_viruses = log_df.shape[1]

    # Fill NaNs with 0 for clustering calculation only
    log_df_cluster = log_df.fillna(0.0)
    if num_viruses > 1:
        linkage_matrix = linkage(log_df_cluster.T, method='average')
        ordered_indices = leaves_list(linkage_matrix)
        # Keep taxonomy pre-sorted order by preserving original column index order if clustering is secondary
        sorted_column_names = [log_df.columns[i] for i in ordered_indices]
        sorted_log_df = log_df[sorted_column_names]
    else:
        sorted_log_df = log_df

    calc_height = max(8.0, min(35.0, num_samples * 0.15 + 4.0))
    calc_width = max(6.0, min(30.0, num_viruses * 0.28 + 3.0))

    plt.figure(figsize=(calc_width, calc_height), dpi=300)

    current_cmap = matplotlib.cm.get_cmap('viridis').copy()
    current_cmap.set_bad(color='#e0e0e0')

    ax = sns.heatmap(sorted_log_df, cmap=current_cmap, linewidths=0.5, linecolor='white', mask=sorted_log_df.isna())

    plt.title(title_text, fontsize=18, pad=15)
    plt.xlabel('Viral Reference / Taxonomy Group', fontsize=14, labelpad=10)
    plt.ylabel('Sample', fontsize=14, labelpad=10)

    xtick_size = max(7, min(12, int(240 / max(1, num_viruses))))
    ytick_size = min(14, max(1, int(300 / max(1, num_samples))))

    plt.xticks(fontsize=xtick_size, rotation=90)
    ax.set_yticks(np.arange(num_samples) + 0.5)
    ax.set_yticklabels(sorted_log_df.index, fontsize=ytick_size, rotation=0)

    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=11)
    for label in colorbar.ax.get_yticklabels():
        label.set_rotation(0)

    plt.savefig(out_path, format="tiff", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[SUCCESS] Heatmap saved to: {out_path}")


def main():
    args = parse_args()

    if not os.path.exists(args.input_report):
        print(f"[ERROR] Master report file not found: {args.input_report}")
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    raw_phase = str(args.phase).strip()
    target_phase = raw_phase if raw_phase and raw_phase.lower() not in ["none", "0"] else "phase"

    raw_project = str(args.project).strip()
    target_project = raw_project if raw_project and raw_project.lower() not in ["none", "0"] else "base"

    target_strategy = args.strategy.strip()

    df = pd.read_csv(args.input_report, sep="\t", dtype=str)

    df.columns = [c.strip() for c in df.columns]
    col_map = {c.lower(): c for c in df.columns}

    sample_col = col_map.get("sample_id", df.columns[0])
    acc_col = col_map.get("virus_accession", df.columns[1])
    reads_col = col_map.get("virus_mapped_reads", df.columns[3])
    cn_col = col_map.get("viral_copy_number", df.columns[8] if len(df.columns) > 8 else df.columns[3])
    name_col = col_map.get("virus_name_sanitized", df.columns[9] if len(df.columns) > 9 else df.columns[1])
    source_col = col_map.get("source_file", df.columns[10] if len(df.columns) > 10 else df.columns[-1])
    phase_col = col_map.get("phase", df.columns[11] if len(df.columns) > 11 else None)
    proj_col = col_map.get("project", df.columns[12] if len(df.columns) > 12 else None)

    sp_taxid_col = col_map.get("species_taxid", "")
    sp_name_col = col_map.get("species_name", "")
    genus_col = col_map.get("genus_name", "")
    family_col = col_map.get("family_name", "")

    # Match phase
    if phase_col:
        df_phase = df[phase_col].astype(str).str.strip()
        mask_phase = (df_phase.str.lower() == target_phase.lower())
        if target_phase.lower() == "phase":
            mask_phase = mask_phase | df_phase.str.lower().isin(["phase", "", "none", "0", "base"])
        df = df[mask_phase]

    # Match project
    if proj_col:
        df_proj = df[proj_col].astype(str).str.strip()
        mask_proj = (df_proj.str.lower() == target_project.lower())
        if target_project.lower() == "base":
            mask_proj = mask_proj | df_proj.str.lower().isin(["base", "", "none", "0"])
        df = df[mask_proj]

    # Match strategy (Source_File)
    if source_col:
        df_source = df[source_col].astype(str).str.strip()
        df_strat = df_source.apply(get_short_strategy)
        mask_strat = (df_strat.str.lower() == target_strategy.lower()) | (df_source.str.lower() == target_strategy.lower())
        if not mask_strat.any():
            mask_strat = df_source.str.lower().str.contains(target_strategy.lower())
        df = df[mask_strat]

    if reads_col in df.columns:
        df_reads_num = pd.to_numeric(df[reads_col], errors='coerce').fillna(0)
        df = df[df_reads_num >= args.min_reads_cutoff].copy()

    if df.empty:
        print(f"[WARNING] No records found with reads >= {args.min_reads_cutoff} for strategy '{target_strategy}'. Skipping heatmap.")
        sys.exit(0)

    value_col = reads_col if args.value_type == "read_counts" else cn_col
    df[value_col] = pd.to_numeric(df[value_col], errors='coerce').fillna(0.0)

    # Build taxonomy attributes for sorting and labeling
    df["tax_family"] = df[family_col].fillna("Unknown").astype(str).str.strip() if family_col else "Unknown"
    df["tax_genus"] = df[genus_col].fillna("Unknown").astype(str).str.strip() if genus_col else "Unknown"
    df["tax_sp_taxid"] = df[sp_taxid_col].fillna("Unknown").astype(str).str.strip() if sp_taxid_col else "Unknown"
    df["tax_sp_name"] = df[sp_name_col].fillna("Unknown").astype(str).str.strip() if sp_name_col else "Unknown"
    df["acc_clean"] = df[acc_col].fillna("Unknown").astype(str).str.strip()
    df["name_clean"] = df[name_col].fillna("Unknown").astype(str).str.strip()

    # Pre-sort df hierarchically by: family_name -> genus_name -> species_taxid -> virus_accession
    df = df.sort_values(by=["tax_family", "tax_genus", "tax_sp_taxid", "acc_clean"])

    # Construct Enhanced Accession-Level Column Display Label:
    # Species_TaxID|Species_Name|Virus_Accession|Virus_Name_Sanitized
    df["Virus_Display"] = (
        df["tax_sp_taxid"] + "|" +
        df["tax_sp_name"] + "|" +
        df["acc_clean"] + "|" +
        df["name_clean"]
    )

    # Preserve sorted unique column order
    unique_cols = list(dict.fromkeys(df["Virus_Display"]))

    pivot_df = df.pivot_table(index=sample_col, columns="Virus_Display", values=value_col, aggfunc="max", fill_value=0.0)
    # Reindex columns to maintain exact hierarchical pre-sorted taxonomic order
    pivot_df = pivot_df.reindex(columns=unique_cols, fill_value=0.0)

    if args.value_type == "copy_number":
        pivot_df = pivot_df[pivot_df.max(axis=1) > 0.0]

    p_val = str(args.phase).strip() if args.phase else ""
    prj_val = str(args.project).strip() if args.project else ""

    is_phase_empty = not p_val or p_val.lower() in ["none", "0", "", "all_cohorts"]
    is_proj_empty = not prj_val or prj_val.lower() in ["none", "0", "", "base", "combined"]

    strat_raw = args.strategy or "clean_flags"
    if "clean_filtered.sorted.flags" in strat_raw or strat_raw == "clean_flags":
        strat_tag = "clean_flags"
    elif "raw_filtered.sorted.flags" in strat_raw or strat_raw == "raw_flags":
        strat_tag = "raw_flags"
    else:
        strat_tag = strat_raw.replace(".cram", "").replace(".bam", "").replace(".", "_")

    if is_phase_empty and is_proj_empty:
        base_out_name = f"{args.dataset}_{args.genome_build}_{strat_tag}_heatmap_{args.value_type}"
    else:
        p_tag = "all" if is_phase_empty else p_val if p_val.startswith("phase") else f"phase{p_val}"
        prj_tag = "base" if is_proj_empty else prj_val
        base_out_name = f"{args.dataset}_{args.genome_build}_{p_tag}_{prj_tag}_{strat_tag}_heatmap_{args.value_type}"

    # 1. Render Accession-Level Heatmap
    out_path_acc = os.path.join(args.out_dir, f"{base_out_name}.tiff")
    title_acc = "Virus Read Counts Heatmap" if args.value_type == "read_counts" else "Virus Copy Number Heatmap"
    render_and_save_heatmap(pivot_df, out_path_acc, title_acc, args.value_type, p_val, prj_val, strat_raw)

    # 2. Grouped Heatmap (if --group-level is set to species, genus, family, or realm)
    grp_level = args.group_level.strip().lower()
    if grp_level != "none":
        col_group_map = {
            "species": df["tax_sp_taxid"] + "|" + df["tax_sp_name"],
            "genus": df["tax_genus"],
            "family": df["tax_family"],
            "realm": df[col_map.get("realm_name", "")].fillna("Unknown").astype(str).str.strip() if "realm_name" in col_map else df["tax_family"]
        }
        group_series = col_group_map.get(grp_level, df["tax_sp_taxid"] + "|" + df["tax_sp_name"])
        df["Group_Display"] = group_series

        unique_groups = list(dict.fromkeys(df["Group_Display"]))
        agg_func = "sum" if args.value_type == "read_counts" else "max"
        
        pivot_grp = df.pivot_table(index=sample_col, columns="Group_Display", values=value_col, aggfunc=agg_func, fill_value=0.0)
        pivot_grp = pivot_grp.reindex(columns=unique_groups, fill_value=0.0)

        if args.value_type == "copy_number":
            pivot_grp = pivot_grp[pivot_grp.max(axis=1) > 0.0]

        out_path_grp = os.path.join(args.out_dir, f"{base_out_name}_grouped_{grp_level}.tiff")
        title_grp = f"Grouped ({grp_level.capitalize()}) Viral {title_acc}"
        render_and_save_heatmap(pivot_grp, out_path_grp, title_grp, args.value_type, p_val, prj_val, strat_raw)


if __name__ == "__main__":
    main()
