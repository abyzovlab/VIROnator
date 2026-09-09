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


def main():
    args = parse_args()

    if not os.path.exists(args.input_report):
        print(f"[ERROR] Master report file not found: {args.input_report}")
        sys.exit(1)

    os.makedirs(args.out_dir, exist_ok=True)

    # 1. Filter Criteria Setup
    # YAML phase: empty -> "phase" in master TSV
    raw_phase = str(args.phase).strip()
    target_phase = raw_phase if raw_phase and raw_phase.lower() not in ["none", "0"] else "phase"

    # YAML project: empty -> "base" in master TSV
    raw_project = str(args.project).strip()
    target_project = raw_project if raw_project and raw_project.lower() not in ["none", "0"] else "base"

    target_strategy = args.strategy.strip()

    # 2. Read Master TSV Report
    # Standard columns:
    # Sample_ID, Virus_Accession, Virus_Length, Virus_Mapped_Reads, Normalized_Coverage,
    # Physical_Coverage, Human_Genome_Size, Sample_Read_Depth, Viral_Copy_Number,
    # Virus_Name_Sanitized, Specimen, phase, project, Source_File
    df = pd.read_csv(args.input_report, sep="\t")

    # Lowercase column names for robust lookup
    df.columns = [c.strip() for c in df.columns]
    col_map = {c.lower(): c for c in df.columns}

    sample_col = col_map.get("sample_id", df.columns[0])
    acc_col = col_map.get("virus_accession", df.columns[1])
    reads_col = col_map.get("virus_mapped_reads", df.columns[3])
    cn_col = col_map.get("viral_copy_number", df.columns[8] if len(df.columns) > 8 else df.columns[3])
    name_col = col_map.get("virus_name_sanitized", df.columns[9] if len(df.columns) > 9 else df.columns[1])
    phase_col = col_map.get("phase", df.columns[11] if len(df.columns) > 11 else None)
    proj_col = col_map.get("project", df.columns[12] if len(df.columns) > 12 else None)
    source_col = col_map.get("source_file", df.columns[13] if len(df.columns) > 13 else df.columns[-1])

    # 3. Apply Filtering
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

    if df.empty:
        print(f"[WARNING] No records found matching phase='{target_phase}', project='{target_project}', strategy='{target_strategy}' in master report.")
        sys.exit(0)

    # 4. Construct Pivot Table
    value_col = reads_col if args.value_type == "read_counts" else cn_col

    # Format Column Names: <Virus_Accession>, <Virus_Name_Sanitized>
    df["Virus_Display"] = df[acc_col].astype(str) + ", " + df[name_col].astype(str)

    # Pivot matrix: Index = Sample_ID, Columns = Virus_Display
    pivot_df = df.pivot_table(index=sample_col, columns="Virus_Display", values=value_col, aggfunc="max", fill_value=0.0)

    if pivot_df.empty or pivot_df.shape[1] == 0:
        print(f"[WARNING] Matrix is empty after pivot for {args.value_type}. Skipping heatmap.")
        sys.exit(0)

    # 5. Transform log10(df + 1e-7) as in OLD/06_heatmap.py
    log_df = np.log10(pivot_df + 1e-7)

    # 6. Hierarchical Clustering (average linkage)
    num_samples = log_df.shape[0]
    num_viruses = log_df.shape[1]

    # Proportional figure dimensions capped for high performance & clean rendering
    calc_height = max(8.0, min(35.0, num_samples * 0.15 + 4.0))
    calc_width = max(10.0, min(30.0, num_viruses * 0.4 + 4.0))

    plt.figure(figsize=(calc_width, calc_height), dpi=300)

    if num_viruses > 1:
        linkage_matrix = linkage(log_df.T, method='average')
        ordered_columns = leaves_list(linkage_matrix)
        sorted_column_names = sorted(log_df.columns[ordered_columns])
        sorted_log_df = log_df[sorted_column_names]
        sorted_log_df = sorted_log_df.loc[:, sorted_column_names]
    else:
        sorted_log_df = log_df

    # 7. Render Heatmap
    ax = sns.heatmap(sorted_log_df, cmap='viridis', linewidths=0.5)

    title_text = 'Virus Read Counts Heatmap' if args.value_type == "read_counts" else 'Virus Copy Number Heatmap'
    plt.title(title_text, fontsize=20, pad=15)
    plt.xlabel('Virus', fontsize=16, labelpad=10)
    plt.ylabel('Sample', fontsize=16, labelpad=10)

    xtick_size = max(6, min(12, int(200 / max(1, num_viruses))))
    ytick_size = max(4, min(10, int(200 / max(1, num_samples))))

    plt.xticks(fontsize=xtick_size, rotation=90)
    plt.yticks(fontsize=ytick_size, rotation=0)

    colorbar = ax.collections[0].colorbar
    colorbar.ax.tick_params(labelsize=12)
    for label in colorbar.ax.get_yticklabels():
        label.set_rotation(0)

    plt.tight_layout()

    # Output filename convention: {dataset}_{genome_build}_heatmap_{value_type}.tiff
    out_filename = f"{args.dataset}_{args.genome_build}_heatmap_{args.value_type}.tiff"
    out_path = os.path.join(args.out_dir, out_filename)
    plt.savefig(out_path, format="tiff", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[SUCCESS] Heatmap successfully generated and saved to: {out_path}")


if __name__ == "__main__":
    main()
