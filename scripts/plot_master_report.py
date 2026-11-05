#!/usr/bin/env python3
"""
scripts/plot_master_report.py
=============================
Standalone, isolated visualization and taxonomy analysis script for VIROnator master reports.

Generates publication-quality figures and statistical tables directly from a master report TSV/CSV:
  1. Taxonomy Statistical Summary Tables (per virus & aggregated by rank: Family, Genus, Species)
  2. Taxonomy Bar Plots across 6 levels (Accession, Species, Genus, Family, Order, Realm)
  3. Family Overlay Plots (Accession bars with shaded family-level background boxes)
  4. Sample Counts & Prevalence Plots (% positive samples & absolute counts, segmented by sample)
  5. Reads per Virus Distribution Plots (Total and sample-level read abundances)
  6. Viral Copy Number Distribution Plots (Total and sample-level copy numbers)
  7. Multi-Virus Density Distributions (KDE density curves of copy numbers for top viruses)
  8. Cross-Sample Viral Abundance Heatmap (Samples vs. Viruses clustered matrix)

Usage:
    python3 scripts/plot_master_report.py --input-report fromNorgan_hg38_master_report.tsv --out-dir ./plots
    python3 scripts/plot_master_report.py --input-report BSMN_hg38_master_report.tsv --strategy clean_flags --format both
"""

import argparse
import os
import pathlib
import re
import sys
import numpy as np
import pandas as pd

# Matplotlib headless backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

# Optional SciPy for Gaussian KDE
try:
    from scipy.stats import gaussian_kde
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# ── Color Palette & Constants ──────────────────────────────────────────────────
BLUE = "#4472C4"

PARENT_GROUP = {
    "Accession": "species_name",
    "Species":   "genus_name",
    "Genus":     "family_name",
    "Family":    "order_name",
    "Order":     "realm_name",
    "Realm":     "realm_name",
}

# Curated 26-color scientific palette (Okabe-Ito + Dark2 core, colorblind-friendly)
PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#7570B3",
    "#56B4E9", "#E7298A", "#66A61E", "#1B9E77", "#CC79A7",
    "#A6761D", "#4B0082", "#FC8D62", "#8DA0CB", "#666666",
    "#E78AC3", "#A6D854", "#984EA3", "#377EB8", "#FFD92F",
    "#E5C494", "#B3DE69", "#BC80BD", "#FB8072", "#80B1D3",
    "#B3B3B3",
]

FAM_PALETTE = ["#C6D8EF", "#FEE8C8", "#E7D4E8", "#FDDBC7", "#E0ECF4", "#E5F5E0", "#FFF7BC", "#D9F0D3"]

BAR_HEIGHT = 0.88
FAM_BOX_MARGIN = 0.46
FONT_TITLE = 11.0
FONT_XLABEL = 10.0
FONT_XTICK = 9.0
FONT_YTICK = 8.5
FONT_FAMILY = 8.5
FONT_LEGEND = 7.5
FONT_LEGEND_TITLE = 8.0
FONT_FOOTNOTE = 8.5


# ── String & Naming Helpers ───────────────────────────────────────────────────
def clean_virus_name(raw_name):
    """Rule-based shortening for standard long viral names."""
    name = str(raw_name).strip()
    name = name.replace("Molluscum contagiosum virus", "MCV").replace("Molluscum_contagiosum_virus", "MCV")
    name = re.sub(r"(?i)human\s+herpesvirus\s+(\d+)", r"HHV-\1", name)
    name = re.sub(r"(?i)torque\s+teno\s+virus\s*(\d*)", r"TTV \1", name)
    name = re.sub(r"(?i)human\s+papillomavirus\s+type\s*(\d+)", r"HPV-\1", name)
    name = re.sub(r"(?i)human\s+papillomavirus\s*(\d+)", r"HPV-\1", name)
    name = re.sub(r"(?i)\s*,?\s*complete\s+genome", "", name)
    return name.strip()


def sanitize_filename(name):
    """Sanitize strings for safe filesystem naming."""
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", str(name)).strip("_")


def detect_axis_break(values, mode):
    """
    Detects whether an extreme outlier gap exists in the positive data distribution.
    Returns (x_cut_low, x_cut_high, x_max_right) if an outlier gap >= 3.5x ratio
    and >= 45% empty gap is detected; otherwise None.
    """
    pos_vals = sorted([float(v) for v in values if v > 0])
    if len(pos_vals) < 2:
        return None

    max_v = pos_vals[-1]
    best_split = None
    for i in range(len(pos_vals) - 1):
        v_low = pos_vals[i]
        v_high = pos_vals[i + 1]
        gap = v_high - v_low
        ratio = v_high / v_low if v_low > 0 else 1.0

        if mode == "reads" and gap < 80:
            continue
        elif mode == "copynumber" and gap < 0.04:
            continue
        elif mode == "samples" and gap < 35.0:
            continue

        if ratio >= 3.5 and (gap / max_v) >= 0.45:
            best_split = (v_low, v_high)
            break

    if not best_split:
        return None

    v_low, v_high = best_split
    if mode == "reads":
        x_cut_low = float(np.ceil(v_low * 1.55 / 10.0) * 10)
        x_cut_high = float(np.floor(v_high * 0.88 / 10.0) * 10)
        x_max_right = max_v * 1.45
    elif mode == "copynumber":
        x_cut_low = v_low * 1.60
        x_cut_high = v_high * 0.85
        x_max_right = max_v * 1.50
    else:  # 'samples' (% positive)
        x_cut_low = float(np.ceil(v_low * 1.50))
        x_cut_high = float(np.floor(v_high * 0.88))
        x_max_right = max_v * 1.45

    if x_cut_low >= x_cut_high:
        return None

    return (x_cut_low, x_cut_high, x_max_right)


def draw_break_markers(ax_left, ax_right, crossing_bars=None, bar_height=0.88):
    """Draws publication-standard chevron notches and zig-zag cut markers on broken axes."""
    fig = ax_left.figure
    fig_w, fig_h = fig.get_size_inches()
    w_l = ax_left.get_position().width * fig_w
    w_r = ax_right.get_position().width * fig_w

    x_l = ax_left.get_position().x1
    x_r = ax_right.get_position().x0
    y_bot = ax_left.get_position().y0
    notch_h_pt = 6.0
    dy = (notch_h_pt / 72.0) / fig_h
    x_mid = (x_l + x_r) / 2.0
    line = mlines.Line2D([x_l, x_mid, x_r], [y_bot, y_bot + dy, y_bot],
                         color="black", linewidth=0.6, transform=fig.transFigure, clip_on=False, zorder=10)
    fig.add_artist(line)

    if crossing_bars:
        amp_pt = 4.2
        dx_l = (amp_pt / 72.0) / w_l
        dx_r = (amp_pt / 72.0) / w_r
        N = 4
        H = bar_height
        for y_bar in crossing_bars:
            y_top = y_bar - H / 2.0
            y_bot_b = y_bar + H / 2.0
            y_pts = [y_top + k * (H / (2 * N)) for k in range(2 * N + 1)]

            verts_l = []
            for k, yp in enumerate(y_pts):
                xp = 1.0 if (k % 2 == 1) else (1.0 - dx_l)
                verts_l.append((xp, yp))
            verts_l.append((1.05, y_bot_b))
            verts_l.append((1.05, y_top))
            poly_l = mpatches.Polygon(verts_l, closed=True, facecolor="white", edgecolor="none",
                                     transform=ax_left.get_yaxis_transform(), clip_on=False, zorder=4)
            ax_left.add_patch(poly_l)

            verts_r = []
            for k, yp in enumerate(y_pts):
                xp = (0.0 + dx_r) if (k % 2 == 1) else 0.0
                verts_r.append((xp, yp))
            verts_r.append((-0.05, y_bot_b))
            verts_r.append((-0.05, y_top))
            poly_r = mpatches.Polygon(verts_r, closed=True, facecolor="white", edgecolor="none",
                                     transform=ax_right.get_yaxis_transform(), clip_on=False, zorder=4)
            ax_right.add_patch(poly_r)


def save_figure(fig, out_dir, filename_base, fmt="both", dpi=300):
    """Saves figure in requested formats (png, tiff, pdf, or both png+tiff)."""
    os.makedirs(out_dir, exist_ok=True)
    formats = []
    if fmt == "both":
        formats = ["png", "tiff"]
    elif fmt == "all":
        formats = ["png", "tiff", "pdf"]
    else:
        formats = [fmt]

    for f in formats:
        out_path = os.path.join(out_dir, f"{filename_base}.{f}")
        try:
            if f == "tiff":
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
            else:
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
            print(f"  ✓ Saved: {out_path}")
        except Exception:
            # Fallback without compression if tiff_lzw is unavailable
            fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
            print(f"  ✓ Saved (standard): {out_path}")


# ── Statistical Table Generation ──────────────────────────────────────────────
def generate_taxonomy_stats_table(df_viral, total_samples, out_dir):
    """
    Computes a comprehensive tabular summary of viral detections across all samples:
    - Virus accession, name, taxonomy ranks
    - Positive samples and % prevalence
    - Reads (Total, Mean, Median, Min, Max)
    - Copy Numbers (Mean, Median, Min, Max)
    - Physical and normalized coverage metrics
    """
    os.makedirs(out_dir, exist_ok=True)
    tsv_out = os.path.join(out_dir, "taxonomy_stats_summary.tsv")
    md_out = os.path.join(out_dir, "taxonomy_stats_summary.md")

    tax_cols = ["virus_accession", "virus_name_sanitized", "species_name", "genus_name", "family_name", "order_name", "realm_name"]
    group_cols = [c for c in tax_cols if c in df_viral.columns]

    records = []
    for keys, sub in df_viral.groupby(group_cols):
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_dict = dict(zip(group_cols, keys))

        n_pos = sub["sample_id"].nunique()
        prev_pct = (n_pos / float(total_samples) * 100.0) if total_samples > 0 else 0.0
        
        reads = sub["virus_mapped_reads"].astype(float)
        cn = sub["viral_copy_number"].astype(float)
        cov = sub["physical_coverage"].astype(float) if "physical_coverage" in sub.columns else pd.Series([0.0])

        rec = {
            "Accession": key_dict.get("virus_accession", "NA"),
            "Virus_Name": key_dict.get("virus_name_sanitized", "NA"),
            "Species": key_dict.get("species_name", "Unknown"),
            "Genus": key_dict.get("genus_name", "Unknown"),
            "Family": key_dict.get("family_name", "Unknown"),
            "Order": key_dict.get("order_name", "Unknown"),
            "Realm": key_dict.get("realm_name", "Unknown"),
            "Positive_Samples": n_pos,
            "Prevalence_Pct": round(prev_pct, 2),
            "Total_Reads": int(reads.sum()),
            "Mean_Reads": round(float(reads.mean()), 1),
            "Median_Reads": round(float(reads.median()), 1),
            "Max_Reads": int(reads.max()),
            "Total_Copy_Number": round(float(cn.sum()), 4),
            "Mean_Copy_Number": round(float(cn.mean()), 4),
            "Median_Copy_Number": round(float(cn.median()), 4),
            "Max_Copy_Number": round(float(cn.max()), 4),
            "Mean_Physical_Cov": round(float(cov.mean()), 2),
        }
        records.append(rec)

    stats_df = pd.DataFrame(records)
    if stats_df.empty:
        return stats_df

    stats_df = stats_df.sort_values(by=["Positive_Samples", "Total_Reads"], ascending=[False, False]).reset_index(drop=True)
    stats_df.to_csv(tsv_out, sep="\t", index=False)
    print(f"  ✓ Written Taxonomy Stats Table: {tsv_out}")

    # Generate Markdown Summary
    with open(md_out, "w") as f_md:
        f_md.write("# VIROnator Viral Taxonomy & Abundance Summary\n\n")
        f_md.write(f"- **Total Cohort Samples**: {total_samples:,}\n")
        f_md.write(f"- **Viral Positive Samples**: {df_viral['sample_id'].nunique():,}\n")
        f_md.write(f"- **Unique Viral Accessions Detected**: {df_viral['virus_accession'].nunique():,}\n")
        f_md.write(f"- **Total Filtered Viral Reads**: {int(df_viral['virus_mapped_reads'].sum()):,}\n\n")
        f_md.write("### Per-Virus Detections (Top 25 by Prevalence)\n\n")
        f_md.write("| Accession | Virus Name | Family | Positive Samples | Prevalence (%) | Total Reads | Mean Reads | Mean Copy Number |\n")
        f_md.write("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")
        for _, r in stats_df.head(25).iterrows():
            f_md.write(f"| {r['Accession']} | {r['Virus_Name']} | {r['Family']} | {r['Positive_Samples']:,} | {r['Prevalence_Pct']:.2f}% | {r['Total_Reads']:,} | {r['Mean_Reads']:.1f} | {r['Mean_Copy_Number']:.4f} |\n")

    print(f"  ✓ Written Markdown Summary: {md_out}")
    return stats_df


# ── Section 1: Taxonomy Level Horizontal Bar Plots ────────────────────────────
def plot_taxonomy_level(df_viral, total_samples, taxid_col, name_col, level_label, out_dir, fmt="both", dpi=300):
    """
    Renders horizontal bar plot for a given taxonomic level (Accession, Species, Genus, Family, Order, Realm)
    colored by parent taxonomic group.
    """
    plot_data = df_viral.copy()
    if level_label != "Accession" and name_col in plot_data.columns:
        plot_data = plot_data[plot_data[name_col] != "Unknown"]

    if plot_data.empty:
        return

    parent_col = PARENT_GROUP[level_label]
    above = []
    if level_label == "Accession":
        above = ["family_name", "genus_name", "species_name"]
    elif level_label == "Species":
        above = ["family_name", "genus_name"]
    elif level_label == "Genus":
        above = ["family_name"]

    group_cols = [c for c in [taxid_col, name_col] if c in plot_data.columns]
    for c in above:
        if c in plot_data.columns and c not in group_cols:
            group_cols.append(c)
    if parent_col in plot_data.columns and parent_col not in group_cols:
        group_cols.append(parent_col)

    grp = (
        plot_data.groupby(group_cols)
        .agg(
            n_positive_samples=("sample_id", "nunique"),
            total_reads=("virus_mapped_reads", "sum"),
        )
        .reset_index()
    )
    if grp.empty:
        return

    grp["pct_positive"] = grp["n_positive_samples"] / float(total_samples) * 100.0 if total_samples > 0 else 0.0
    label_col = name_col if name_col in grp.columns else group_cols[0]
    if taxid_col in grp.columns and taxid_col != name_col:
        grp["label"] = grp[taxid_col].astype(str) + ": " + grp[name_col].astype(str)
    else:
        grp["label"] = grp[label_col].astype(str)

    # Hierarchical sorting
    sort_cols, asc_flags = [], []
    for c in above:
        if c in plot_data.columns:
            c_tot = plot_data.groupby(c)["sample_id"].nunique().to_dict()
            grp[c + "_tot"] = grp[c].map(c_tot)
            grp[c + "_unk"] = grp[c] == "Unknown"
            sort_cols.extend([c + "_unk", c + "_tot", c])
            asc_flags.extend([True, False, True])

    sort_cols.extend(["n_positive_samples", "total_reads"])
    asc_flags.extend([False, False])
    grp = grp.sort_values(by=sort_cols, ascending=asc_flags).reset_index(drop=True)

    # Assign parent colors
    if parent_col in grp.columns:
        unique_parents = list(dict.fromkeys(grp[parent_col]))
        color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(unique_parents)}
        bar_colors = [color_map[g] for g in grp[parent_col]]
    else:
        unique_parents = ["Virus"]
        color_map = {"Virus": BLUE}
        bar_colors = [BLUE] * len(grp)

    # Plotting
    n_bars = len(grp)
    fig_h = max(3.5, n_bars * 0.38)
    fig, ax = plt.subplots(figsize=(10.5, fig_h))

    bars = ax.barh(grp["label"], grp["pct_positive"], color=bar_colors, edgecolor="none", height=BAR_HEIGHT)
    ax.invert_yaxis()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    max_val = grp["pct_positive"].max() if len(grp) else 1
    for bar, reads in zip(bars, grp["total_reads"]):
        w = bar.get_width()
        if w > max_val * 0.30:
            ax.text(w - max_val * 0.015, bar.get_y() + bar.get_height() / 2, f"{int(reads):,} reads",
                    va="center", ha="right", fontsize=FONT_YTICK, color="white")
        else:
            ax.text(w + max_val * 0.015, bar.get_y() + bar.get_height() / 2, f"{int(reads):,} reads",
                    va="center", ha="left", fontsize=FONT_YTICK, color="black")

    if parent_col in grp.columns and len(unique_parents) > 1:
        parent_title = parent_col.replace("_name", "").replace("_", " ").title()
        handles = [mpatches.Patch(facecolor=color_map[g], edgecolor="none", label=g) for g in unique_parents]
        ax.legend(handles=handles, title=parent_title, fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE,
                  bbox_to_anchor=(1.01, 0.5), loc="center left", frameon=False)

    ax.set_xlabel(f"% positive samples (n = {total_samples:,})", fontsize=FONT_XLABEL)
    ax.set_title(f"Distribution of Viruses – {level_label} Level", fontsize=FONT_TITLE)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
    ax.set_xlim(0, max_val * 1.35)
    ax.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
    ax.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

    plt.tight_layout()
    save_figure(fig, out_dir, f"taxonomy_{level_label.lower()}", fmt=fmt, dpi=dpi)
    plt.close(fig)


# ── Section 2: Family Overlay Horizontal Bar Plot ─────────────────────────────
def plot_family_overlay(df_viral, total_samples, out_dir, fmt="both", dpi=300):
    """
    Renders Accession-level horizontal bars with translucent shaded family background blocks.
    """
    if "family_name" not in df_viral.columns:
        return

    grp = (
        df_viral.groupby(["virus_accession", "virus_name_sanitized", "family_name", "genus_name", "species_name"])
        .agg(
            n_positive_samples=("sample_id", "nunique"),
            total_reads=("virus_mapped_reads", "sum"),
        )
        .reset_index()
    )
    if grp.empty:
        return

    grp["pct_positive"] = grp["n_positive_samples"] / float(total_samples) * 100.0 if total_samples > 0 else 0.0
    grp["label"] = grp["virus_accession"].astype(str) + ": " + grp["virus_name_sanitized"].astype(str)

    # Sort strictly by Family total samples -> Genus -> Species -> Accession
    fam_tot = df_viral.groupby("family_name")["sample_id"].nunique().to_dict()
    grp["fam_tot"] = grp["family_name"].map(fam_tot)
    grp["fam_unk"] = grp["family_name"] == "Unknown"

    grp = grp.sort_values(by=["fam_unk", "fam_tot", "family_name", "n_positive_samples", "total_reads"],
                          ascending=[True, False, True, False, False]).reset_index(drop=True)

    n_bars = len(grp)
    fig_h = max(4.0, n_bars * 0.40)
    fig, ax = plt.subplots(figsize=(11.0, fig_h))

    # Family background shaded blocks
    fam_blocks = []
    current_fam = None
    start_idx = 0
    for idx, fam in enumerate(grp["family_name"]):
        if fam != current_fam:
            if current_fam is not None:
                fam_blocks.append((current_fam, start_idx, idx - 1))
            current_fam = fam
            start_idx = idx
    if current_fam is not None:
        fam_blocks.append((current_fam, start_idx, len(grp) - 1))

    max_val = grp["pct_positive"].max()
    fam_pct_tot = {f: (fam_tot[f] / float(total_samples) * 100.0) for f in fam_tot}
    global_max_x = max(max_val, max(fam_pct_tot.values()) if fam_pct_tot else max_val)

    for i, (fam, s_idx, e_idx) in enumerate(fam_blocks):
        tint = FAM_PALETTE[i % len(FAM_PALETTE)]
        w_fam = fam_pct_tot.get(fam, 0.0)
        y_top = s_idx - FAM_BOX_MARGIN
        box_h = (e_idx - s_idx + 1) - (1.0 - 2 * FAM_BOX_MARGIN)
        rect = mpatches.Rectangle((0, y_top), w_fam, box_h, facecolor=tint, edgecolor="none", alpha=0.55, zorder=1)
        ax.add_patch(rect)

        # Family label
        mid_y = (s_idx + e_idx) / 2.0
        ax.text(w_fam + global_max_x * 0.015, mid_y, f"{fam} ({w_fam:.1f}%)",
                va="center", ha="left", fontsize=FONT_FAMILY, fontweight="bold", color="#333333", zorder=3)

    bars = ax.barh(grp["label"], grp["pct_positive"], color=BLUE, edgecolor="none", height=BAR_HEIGHT, zorder=2)
    ax.invert_yaxis()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    for bar, reads in zip(bars, grp["total_reads"]):
        w = bar.get_width()
        if w > global_max_x * 0.35:
            ax.text(w - global_max_x * 0.015, bar.get_y() + bar.get_height() / 2, f"{int(reads):,} reads",
                    va="center", ha="right", fontsize=FONT_YTICK, color="white", zorder=4)
        else:
            ax.text(w + global_max_x * 0.015, bar.get_y() + bar.get_height() / 2, f"{int(reads):,} reads",
                    va="center", ha="left", fontsize=FONT_YTICK, color="black", zorder=4)

    ax.set_xlabel(f"% positive samples (n = {total_samples:,})", fontsize=FONT_XLABEL)
    ax.set_title("Distribution of Viruses – Accession Level (Family Overlay)", fontsize=FONT_TITLE)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
    ax.set_xlim(0, global_max_x * 1.45)
    ax.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
    ax.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

    plt.tight_layout()
    save_figure(fig, out_dir, "taxonomy_accession_family_overlay", fmt=fmt, dpi=dpi)
    plt.close(fig)


# ── Section 3: Segmented Sample-Level Enrichment (Counts, Reads, Copy Number) ─
def plot_sample_enrichment_segmented(df_viral, total_samples, mode, out_dir, fmt="both", dpi=300):
    """
    Renders stacked/segmented horizontal bars per virus colored by individual sample IDs:
    - mode="samples"    : % positive samples and sample count
    - mode="reads"      : total and per-sample mapped reads
    - mode="copynumber" : total and per-sample viral copy numbers
    """
    if "family_name" not in df_viral.columns:
        return

    grp = (
        df_viral.groupby(["virus_accession", "virus_name_sanitized", "family_name"])
        .agg(
            n_positive_samples=("sample_id", "nunique"),
            total_reads=("virus_mapped_reads", "sum"),
            total_copynumber=("viral_copy_number", "sum"),
        )
        .reset_index()
    )
    if grp.empty:
        return

    grp["pct_positive"] = grp["n_positive_samples"] / float(total_samples) * 100.0 if total_samples > 0 else 0.0
    grp["label"] = grp["virus_accession"].astype(str) + ": " + grp["virus_name_sanitized"].astype(str)

    # Sort
    fam_tot = df_viral.groupby("family_name")["sample_id"].nunique().to_dict()
    grp["fam_tot"] = grp["family_name"].map(fam_tot)
    if mode == "samples":
        grp = grp.sort_values(by=["fam_tot", "family_name", "n_positive_samples", "total_reads"],
                              ascending=[False, True, False, False]).reset_index(drop=True)
    elif mode == "copynumber":
        grp = grp.sort_values(by=["fam_tot", "family_name", "total_copynumber", "n_positive_samples"],
                              ascending=[False, True, False, False]).reset_index(drop=True)
    else:
        grp = grp.sort_values(by=["fam_tot", "family_name", "total_reads", "n_positive_samples"],
                              ascending=[False, True, False, False]).reset_index(drop=True)

    unique_samples = list(dict.fromkeys(df_viral["sample_id"]))
    sample_color_map = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(unique_samples)}

    n_bars = len(grp)
    fig_h = max(4.5, n_bars * 0.42)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))

    # Family background tinting
    fam_order = list(dict.fromkeys(grp["family_name"]))
    fam_tint_map = {f: FAM_PALETTE[i % len(FAM_PALETTE)] for i, f in enumerate(fam_order)}

    # Segmented bars
    for idx, row in grp.iterrows():
        acc = row["virus_accession"]
        sub_samples = df_viral[df_viral["virus_accession"] == acc]

        cur_left = 0.0
        for _, s_row in sub_samples.iterrows():
            s_id = s_row["sample_id"]
            if mode == "samples":
                seg_w = 1.0 / float(total_samples) * 100.0 if total_samples > 0 else 1.0
            elif mode == "copynumber":
                seg_w = float(s_row["viral_copy_number"])
            else:
                seg_w = float(s_row["virus_mapped_reads"])

            c = sample_color_map[s_id]
            ax.barh(idx, seg_w, left=cur_left, color=c, edgecolor="none", height=BAR_HEIGHT, zorder=3)
            cur_left += seg_w

        # Total annotation text
        if mode == "samples":
            txt = f"{row['n_positive_samples']:,} samples ({row['pct_positive']:.1f}%)"
        elif mode == "copynumber":
            txt = f"CN: {row['total_copynumber']:.3f} ({row['n_positive_samples']} samples)"
        else:
            txt = f"{int(row['total_reads']):,} reads ({row['n_positive_samples']} samples)"

        ax.text(cur_left * 1.015 + 0.1, idx, txt, va="center", ha="left", fontsize=FONT_YTICK, color="black", zorder=4)

    ax.set_yticks(range(len(grp)))
    ax.set_yticklabels(grp["label"], fontsize=FONT_YTICK)
    ax.invert_yaxis()

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    if mode == "samples":
        ax.set_xlabel(f"% positive samples (n = {total_samples:,})", fontsize=FONT_XLABEL)
        ax.set_title("Viral Sample Counts & Prevalence (Segmented by Sample)", fontsize=FONT_TITLE)
        ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
        fname = "viral_sample_prevalence"
    elif mode == "copynumber":
        ax.set_xlabel("Total Viral Copy Number across Cohort", fontsize=FONT_XLABEL)
        ax.set_title("Viral Copy Number Distribution (Segmented by Sample)", fontsize=FONT_TITLE)
        fname = "viral_copy_number_distribution"
    else:
        ax.set_xlabel("Total Mapped Viral Reads", fontsize=FONT_XLABEL)
        ax.set_title("Viral Mapped Reads Distribution (Segmented by Sample)", fontsize=FONT_TITLE)
        fname = "viral_mapped_reads_distribution"

    ax.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
    ax.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

    plt.tight_layout()
    save_figure(fig, out_dir, fname, fmt=fmt, dpi=dpi)
    plt.close(fig)


# ── Section 4: Multi-Virus Copy Number Density Distributions (KDE) ────────────
def plot_multi_virus_density(df_viral, out_dir, top_n=10, fmt="both", dpi=300):
    """
    Renders smoothed continuous probability density curves (KDE) of viral copy numbers
    for the top N most prevalent viruses.
    """
    if "viral_copy_number" not in df_viral.columns or df_viral.empty:
        return

    top_viruses = (
        df_viral.groupby("virus_name_sanitized")["sample_id"]
        .nunique()
        .sort_values(ascending=False)
        .head(top_n)
        .index.tolist()
    )

    fig, ax = plt.subplots(figsize=(10.0, 5.5))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    plotted_any = False
    for i, virus in enumerate(top_viruses):
        cn_vals = df_viral[df_viral["virus_name_sanitized"] == virus]["viral_copy_number"].astype(float).values
        cn_vals = cn_vals[cn_vals > 0]
        if len(cn_vals) < 2:
            continue

        color = PALETTE[i % len(PALETTE)]
        label_str = f"{virus} (n = {len(cn_vals)})"

        if HAS_SCIPY:
            try:
                log_vals = np.log10(cn_vals)
                kde = gaussian_kde(log_vals)
                x_eval = np.linspace(min(log_vals) - 0.5, max(log_vals) + 0.5, 300)
                y_eval = kde(x_eval)
                ax.plot(x_eval, y_eval, label=label_str, color=color, linewidth=2.0)
                ax.fill_between(x_eval, y_eval, alpha=0.12, color=color)
                plotted_any = True
            except Exception:
                pass
        else:
            # Fallback to histogram if scipy is not installed
            counts, bin_edges = np.histogram(np.log10(cn_vals), bins=15, density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            ax.plot(bin_centers, counts, label=label_str, color=color, linewidth=2.0)
            plotted_any = True

    if not plotted_any:
        plt.close(fig)
        return

    ax.set_xlabel(r"$\log_{10}$ (Viral Copy Number)", fontsize=FONT_XLABEL)
    ax.set_ylabel("Probability Density", fontsize=FONT_XLABEL)
    ax.set_title(f"Multi-Virus Copy Number Density Distribution (Top {len(top_viruses)} Viruses)", fontsize=FONT_TITLE)
    ax.legend(fontsize=FONT_LEGEND, title="Virus", title_fontsize=FONT_LEGEND_TITLE, frameon=False)
    ax.tick_params(axis="both", labelsize=FONT_XTICK)

    plt.tight_layout()
    save_figure(fig, out_dir, "viral_copy_number_multi_density", fmt=fmt, dpi=dpi)
    plt.close(fig)


# ── Section 5: Cross-Sample Viral Abundance Heatmap ───────────────────────────
def plot_abundance_heatmap(df_viral, out_dir, top_n_viruses=25, fmt="both", dpi=300):
    """
    Renders cross-sample clustered heatmap of viral copy numbers across top viruses.
    """
    if "viral_copy_number" not in df_viral.columns or df_viral.empty:
        return

    top_v = (
        df_viral.groupby("virus_name_sanitized")["sample_id"]
        .nunique()
        .sort_values(ascending=False)
        .head(top_n_viruses)
        .index.tolist()
    )
    df_sub = df_viral[df_viral["virus_name_sanitized"].isin(top_v)]

    # Pivot matrix: Viruses (rows) x Samples (columns)
    matrix = df_sub.pivot_table(index="virus_name_sanitized", columns="sample_id",
                                values="viral_copy_number", fill_value=0.0)

    if matrix.empty:
        return

    # Log10 transform with pseudocount
    log_matrix = np.log10(matrix + 1e-4)

    fig, ax = plt.subplots(figsize=(max(8.0, len(matrix.columns) * 0.12), max(5.0, len(matrix) * 0.35)))
    cax = ax.imshow(log_matrix.values, cmap="YlOrRd", aspect="auto", interpolation="nearest")

    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index, fontsize=FONT_YTICK)

    # Show x-ticks only if <= 40 samples to prevent crowding
    if len(matrix.columns) <= 40:
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_xticklabels(matrix.columns, rotation=90, fontsize=FONT_XTICK)
    else:
        ax.set_xticks([])
        ax.set_xlabel(f"Cohort Samples (n = {len(matrix.columns):,})", fontsize=FONT_XLABEL)

    ax.set_title(f"Viral Copy Number Heatmap (Top {len(matrix.index)} Viruses)", fontsize=FONT_TITLE)
    cbar = fig.colorbar(cax, ax=ax, fraction=0.03, pad=0.03)
    cbar.set_label(r"$\log_{10}(\mathrm{Copy\ Number} + 10^{-4})$", fontsize=FONT_XTICK)

    plt.tight_layout()
    save_figure(fig, out_dir, "viral_abundance_heatmap", fmt=fmt, dpi=dpi)
    plt.close(fig)


# ── Main Controller ───────────────────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(
        description="Isolated, standalone script to generate taxonomy statistics and publication plots from a VIROnator master report."
    )
    parser.add_argument("--input-report", required=True,
                        help="Path to consolidated master report TSV/CSV file (e.g. fromNorgan_hg38_master_report.tsv)")
    parser.add_argument("--out-dir", default="./master_report_plots",
                        help="Output directory to save plots and summary tables (default: ./master_report_plots)")
    parser.add_argument("--strategy", default=None,
                        help="Optional filter for source_file CRAM strategy (e.g. exogeneSR_viral_clean_filtered.sorted.flags.cram)")
    parser.add_argument("--phase", default=None,
                        help="Optional filter for specific cohort phase (e.g. panviral, hpv, phase1)")
    parser.add_argument("--project", default=None,
                        help="Optional filter for specific project (e.g. base)")
    parser.add_argument("--short-names", default=None,
                        help="Optional path to short names TSV mapping file (final_clean_id \\t short_name)")
    parser.add_argument("--format", default="both", choices=["png", "tiff", "pdf", "both", "all"],
                        help="Image format for figures: 'png', 'tiff', 'pdf', 'both' (png+tiff), or 'all' (default: both)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Plot resolution DPI (default: 300)")
    parser.add_argument("--top-n", type=int, default=15,
                        help="Top N viruses for density and heatmap plots (default: 15)")
    return parser.parse_args()


def main():
    args = parse_args()

    report_path = pathlib.Path(args.input_report)
    if not report_path.exists():
        print(f"[FATAL ERROR] Master report file not found at: {report_path}", file=sys.stderr)
        sys.exit(1)

    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"==================================================================")
    print(f" Loading VIROnator Master Report: {report_path.name}")
    print(f" Destination Directory:          {out_dir}")
    print(f"==================================================================")

    # Detect delimiter
    sep = "\t" if report_path.suffix.lower() in [".tsv", ".txt"] else ","
    df = pd.read_csv(report_path, sep=sep, low_memory=False)
    print(f"Total records in report  : {len(df):,}")

    # Standardize column types
    if "viral_copy_number" in df.columns:
        df["viral_copy_number"] = pd.to_numeric(df["viral_copy_number"], errors="coerce").fillna(0.0)
    if "virus_mapped_reads" in df.columns:
        df["virus_mapped_reads"] = pd.to_numeric(df["virus_mapped_reads"], errors="coerce").fillna(0)

    # Apply strategy filter if specified
    if args.strategy:
        mask = df["source_file"].str.contains(args.strategy, case=False, na=False)
        df_filt = df[mask].copy()
        print(f"Filtered by strategy '{args.strategy}': {len(df_filt):,} rows remaining")
    else:
        # Default priority: clean_filtered flags CRAM if present, otherwise all
        clean_mask = df["source_file"].str.contains("clean_filtered.*flags", case=False, na=False)
        if clean_mask.any():
            df_filt = df[clean_mask].copy()
            print(f"Auto-selected primary clean filtered strategy: {len(df_filt):,} rows")
        else:
            df_filt = df.copy()

    # Apply phase and project filters if specified
    if args.phase:
        p_clean = args.phase.replace("phase", "").lower()
        df_filt = df_filt[df_filt["phase"].astype(str).str.lower().str.replace("phase", "") == p_clean].copy()
        print(f"Filtered by phase '{args.phase}': {len(df_filt):,} rows remaining")
    if args.project:
        df_filt = df_filt[df_filt["project"].astype(str).str.lower() == args.project.lower()].copy()
        print(f"Filtered by project '{args.project}': {len(df_filt):,} rows remaining")

    total_cohort_samples = df_filt["sample_id"].nunique()
    print(f"Total unique cohort samples in subset: {total_cohort_samples:,}")

    # Exclude non-viral / empty accession rows for viral analyses
    viral_mask = (
        df_filt["virus_accession"].notna() &
        (df_filt["virus_accession"] != "") &
        (df_filt["virus_accession"].astype(str).str.lower() != "none") &
        (df_filt["virus_mapped_reads"] > 0)
    )
    df_viral = df_filt[viral_mask].copy()
    print(f"Positive viral detection rows         : {len(df_viral):,}")
    print(f"Viral positive samples                : {df_viral['sample_id'].nunique():,}")
    print(f"Unique viral taxa                     : {df_viral['virus_accession'].nunique():,}")

    if df_viral.empty:
        print("[WARNING] No viral positive records found after applying filters. Exiting.")
        sys.exit(0)

    # Optional short names mapping
    if args.short_names and os.path.exists(args.short_names):
        try:
            sm_df = pd.read_csv(args.short_names, sep="\t")
            id_col = sm_df.columns[0]
            name_col = sm_df.columns[1] if len(sm_df.columns) > 1 else id_col
            s_map = dict(zip(sm_df[id_col].astype(str).str.strip(), sm_df[name_col].astype(str).str.strip()))
            df_viral["virus_name_sanitized"] = df_viral["virus_accession"].map(s_map).fillna(df_viral["virus_name_sanitized"])
            print(f"Loaded short names map ({len(s_map):,} entries)")
        except Exception as e:
            print(f"[WARNING] Could not parse short names file: {e}")

    # Standardize virus names
    if "virus_name_sanitized" in df_viral.columns:
        df_viral["virus_name_sanitized"] = df_viral["virus_name_sanitized"].apply(clean_virus_name)

    # 1. Generate Comprehensive Statistical Tables
    print("\n── 1. Generating Taxonomy Statistical Summary Tables ──")
    generate_taxonomy_stats_table(df_viral, total_cohort_samples, out_dir)

    # 2. Generate Taxonomy Horizontal Bar Plots (6 Levels)
    print("\n── 2. Generating Taxonomy Level Bar Plots ──")
    tax_levels = [
        ("virus_accession", "virus_name_sanitized", "Accession"),
        ("species_taxid", "species_name", "Species"),
        ("genus_taxid", "genus_name", "Genus"),
        ("family_taxid", "family_name", "Family"),
        ("order_taxid", "order_name", "Order"),
        ("realm_taxid", "realm_name", "Realm"),
    ]
    for taxid_col, name_col, lvl_lbl in tax_levels:
        plot_taxonomy_level(df_viral, total_cohort_samples, taxid_col, name_col, lvl_lbl, out_dir, fmt=args.format, dpi=args.dpi)

    # 3. Generate Family Background Overlay Plot
    print("\n── 3. Generating Accession Family Overlay Plot ──")
    plot_family_overlay(df_viral, total_cohort_samples, out_dir, fmt=args.format, dpi=args.dpi)

    # 4. Generate Segmented Sample Enrichment Plots (Prevalence, Reads, Copy Number)
    print("\n── 4. Generating Sample Counts, Reads, and Copy Number Plots ──")
    plot_sample_enrichment_segmented(df_viral, total_cohort_samples, mode="samples", out_dir=out_dir, fmt=args.format, dpi=args.dpi)
    plot_sample_enrichment_segmented(df_viral, total_cohort_samples, mode="reads", out_dir=out_dir, fmt=args.format, dpi=args.dpi)
    plot_sample_enrichment_segmented(df_viral, total_cohort_samples, mode="copynumber", out_dir=out_dir, fmt=args.format, dpi=args.dpi)

    # 5. Generate Multi-Virus Copy Number Density KDE Plot
    print("\n── 5. Generating Multi-Virus Copy Number Density Distribution ──")
    plot_multi_virus_density(df_viral, out_dir=out_dir, top_n=args.top_n, fmt=args.format, dpi=args.dpi)

    # 6. Generate Viral Abundance Heatmap
    print("\n── 6. Generating Cross-Sample Abundance Heatmap ──")
    plot_abundance_heatmap(df_viral, out_dir=out_dir, top_n_viruses=25, fmt=args.format, dpi=args.dpi)

    print(f"\n==================================================================")
    print(f" [DONE] All plots and statistical summaries successfully saved to:")
    print(f" {out_dir.resolve()}")
    print(f"==================================================================")


if __name__ == "__main__":
    main()
