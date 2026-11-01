#!/usr/bin/env python3
"""
scripts/generate_stats.py
=========================
Unified Stats & Visualizations Module (Module 5)
Integrates:
  5.1 Core Summary Tables & Reports
  5.2 Taxonomic Distribution & Landscape Plots (Broken & Non-Broken Axes, Family Background Overlays)
  5.3 Cohort Demographics Plots (Prioritized from OLD/taxonomy_demographics_plots.py, Non-Broken Axes)
  5.4 Abundance & Multi-Virus Density Distributions (Merged from Module 7, Smoothed Multi-Virus Curves)
  5.5 Cross-Sample Heatmaps & Ubiquitous Virus Reporting
  5.6 Spatial Uniformity & Entropy Plots (H and J comparative multi-virus boxplots)
  5.7 Coverage Depth, Hotspots & Artifact Exclusion (Pooled & Per-sample in coverage_plots/)
"""

import argparse
import collections
import glob
import os
import pathlib
import re
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

# ── Palette & Typography Constants from OLD/taxonomy_demographics_plots.py ──────
BLUE = "#4472C4"

PARENT_GROUP = {
    "Accession": "species_name",
    "Species":   "genus_name",
    "Genus":     "family_name",
    "Family":    "order_name",
    "Order":     "realm_name",
    "Realm":     "realm_name",
}

PALETTE = [
    "#0072B2", "#E69F00", "#009E73", "#D55E00", "#7570B3",
    "#56B4E9", "#E7298A", "#66A61E", "#1B9E77", "#CC79A7",
    "#A6761D", "#4B0082", "#FC8D62", "#8DA0CB", "#666666",
    "#E78AC3", "#A6D854", "#984EA3", "#377EB8", "#FFD92F",
    "#E5C494", "#B3DE69", "#BC80BD", "#FB8072", "#80B1D3",
    "#B3B3B3",
]

HIERARCHY = ["realm_name", "order_name", "family_name", "genus_name", "species_name"]

BAR_HEIGHT = 0.88
FAM_BOX_MARGIN = 0.46
FONT_TITLE = 10.5
FONT_XLABEL = 10.0
FONT_XTICK = 9.0
FONT_YTICK = 8.5
FONT_FAMILY = 8.5
FONT_LEGEND = 7.5
FONT_LEGEND_TITLE = 8.0
FONT_FOOTNOTE = 9.0

CLASS_COLORS = {
    "Saliva": "#1f77b4",
    "Blood":  "#d62728",
    "Brain":  "#2ca02c",
    "Other":  "#7f7f7f",
}
CLASS_ORDER = ["Saliva", "Blood", "Brain", "Other"]


def clean_virus_name(raw_name):
    """Fallback automated rule-based shortening."""
    name = str(raw_name).strip()
    name = name.replace("Molluscum contagiosum virus", "MCV").replace("Molluscum_contagiosum_virus", "MCV")
    name = re.sub(r"(?i)human\s+herpesvirus\s+(\d+)", r"HHV-", name)
    name = re.sub(r"(?i)torque\s+teno\s+virus\s*(\d*)", r"TTV", name)
    name = re.sub(r"(?i)human\s+papillomavirus\s+type\s*(\d+)", r"HPV-", name)
    name = re.sub(r"(?i)human\s+papillomavirus\s*(\d+)", r"HPV-", name)
    name = re.sub(r"(?i)\s*,?\s*complete\s+genome", "", name)
    return name


def sanitize_filename(name):
    return re.sub(r"[^A-Za-z0-9_\-]+", "_", str(name)).strip("_")


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
    return sanitize_filename(fname)


def detect_axis_break(values, mode):
    """
    Detects whether an extreme outlier gap exists in the positive data distribution.
    If an outlier with >= 3.5x ratio and >= 45% empty gap of the axis range is found,
    returns (x_cut_low, x_cut_high, x_max_right).
    Otherwise, returns None (continuous linear scale, no break).
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

        # Minimum absolute gap thresholds to avoid splitting on tiny differences
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
    """Draws publication-standard zig-zag break markers across the broken bottom axis spine and across any horizontal bars crossing the cut, with no additional lines."""
    fig = ax_left.figure
    fig_w, fig_h = fig.get_size_inches()
    w_l = ax_left.get_position().width * fig_w
    w_r = ax_right.get_position().width * fig_w

    # 1. Chevron notch on bottom axis spine (pointing upward towards the bars)
    x_l = ax_left.get_position().x1
    x_r = ax_right.get_position().x0
    y_bot = ax_left.get_position().y0
    notch_h_pt = 6.0
    dy = (notch_h_pt / 72.0) / fig_h
    x_mid = (x_l + x_r) / 2.0
    line = mlines.Line2D([x_l, x_mid, x_r], [y_bot, y_bot + dy, y_bot], color="black", linewidth=0.6, transform=fig.transFigure, clip_on=False, zorder=10)
    fig.add_artist(line)

    # 2. Zig-Zag cut for crossing bars (no extra lines, pure zig-zag bar ends)
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

            # Left mask (on ax_left): covers to the right of the zig-zag profile
            verts_l = []
            for k, yp in enumerate(y_pts):
                xp = 1.0 if (k % 2 == 1) else (1.0 - dx_l)
                verts_l.append((xp, yp))
            verts_l.append((1.05, y_bot_b))
            verts_l.append((1.05, y_top))
            poly_l = mpatches.Polygon(verts_l, closed=True, facecolor="white", edgecolor="none", transform=ax_left.get_yaxis_transform(), clip_on=False, zorder=4)
            ax_left.add_patch(poly_l)

            # Right mask (on ax_right): covers to the left of the matching zig-zag profile
            verts_r = []
            for k, yp in enumerate(y_pts):
                xp = (0.0 + dx_r) if (k % 2 == 1) else 0.0
                verts_r.append((xp, yp))
            verts_r.append((-0.05, y_bot_b))
            verts_r.append((-0.05, y_top))
            poly_r = mpatches.Polygon(verts_r, closed=True, facecolor="white", edgecolor="none", transform=ax_right.get_yaxis_transform(), clip_on=False, zorder=4)
            ax_right.add_patch(poly_r)



def plot_taxonomy_level(data, taxid_col, name_col, level_label, save_path_prefix, n_total_cohort=0, strat_tag="", enable_broken=True):
    """
    Plots horizontal bars for virus taxonomy level colored by parent group.
    Generates both nonbroken and (if detected) broken-axis plots.
    """
    plot_data = data.copy()
    if level_label != "Accession":
        plot_data = plot_data[plot_data[name_col] != "Unknown"]

    parent_col = PARENT_GROUP[level_label]
    if level_label == "Accession":
        above = ["family_name", "genus_name", "species_name"]
    elif level_label == "Species":
        above = ["family_name", "genus_name"]
    elif level_label == "Genus":
        above = ["family_name"]
    else:
        above = []

    group_cols = [taxid_col, name_col]
    for c in above:
        if c not in group_cols:
            group_cols.append(c)
    if parent_col not in group_cols and parent_col in plot_data.columns:
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

    n_tot = n_total_cohort if n_total_cohort > 0 else plot_data["sample_id"].nunique()
    grp["pct_positive"] = grp["n_positive_samples"] / float(n_tot) * 100.0 if n_tot > 0 else 0.0
    grp["label"] = grp[taxid_col].astype(str) + ": " + grp[name_col].astype(str)

    sort_cols = []
    asc_flags = []
    for c in above:
        c_tot = plot_data.groupby(c)["sample_id"].nunique().to_dict()
        grp[c + "_tot"] = grp[c].map(c_tot)
        grp[c + "_unk"] = grp[c] == "Unknown"
        sort_cols.extend([c + "_unk", c + "_tot", c])
        asc_flags.extend([True, False, True])

    sort_cols.extend(["n_positive_samples", "total_reads"])
    asc_flags.extend([False, False])
    grp = grp.sort_values(by=sort_cols, ascending=asc_flags).reset_index(drop=True)

    unique_groups = list(dict.fromkeys(grp[parent_col]))
    color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(unique_groups)}
    bar_colors = [color_map[g] for g in grp[parent_col]]

    axis_break_detected = detect_axis_break(grp["pct_positive"], "samples") if enable_broken else None
    break_variants = []
    if axis_break_detected:
        break_variants.append((axis_break_detected, f"{save_path_prefix}__brokenaxis.tiff"))
    break_variants.append((None, f"{save_path_prefix}__nonbrokenaxis.tiff"))

    n_bars = len(grp)
    fig_h = max(3.5, n_bars * 0.40)

    for axis_break, save_path in break_variants:
        if axis_break:
            x_cut_low, x_cut_high, x_max_right = axis_break
            fig, (ax, ax_right) = plt.subplots(1, 2, figsize=(11.5, fig_h), sharey=True, gridspec_kw={"width_ratios": [3.2, 1.0], "wspace": 0.045})
            axes_list = [ax, ax_right]
            for a in axes_list:
                a.spines["top"].set_visible(False)
                a.spines["right"].set_visible(False)
                a.spines["left"].set_linewidth(0.6)
                a.spines["bottom"].set_linewidth(0.6)
            ax_right.spines["left"].set_visible(False)
            ax_right.tick_params(left=False)
        else:
            fig, ax = plt.subplots(figsize=(10.5, fig_h))
            axes_list = [ax]
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(0.6)
            ax.spines["bottom"].set_linewidth(0.6)

        crossing_bars = []
        for a in axes_list:
            bars = a.barh(grp["label"], grp["pct_positive"], color=bar_colors, edgecolor="none", height=BAR_HEIGHT)
        ax.invert_yaxis()

        max_val = grp["pct_positive"].max() if len(grp) else 1
        for idx, (val, reads) in enumerate(zip(grp["pct_positive"], grp["total_reads"])):
            if axis_break and val > x_cut_low:
                crossing_bars.append(idx)
                target_ax = ax_right
                x_text = val + (x_max_right - x_cut_high) * 0.02
                target_ax.text(x_text, idx, f"{int(reads):,} reads", va="center", ha="left", fontsize=FONT_YTICK, color="black", zorder=5)
            else:
                target_ax = ax
                if val > max_val * 0.25 and not axis_break:
                    target_ax.text(val - max_val * 0.015, idx, f"{int(reads):,} reads", va="center", ha="right", fontsize=FONT_YTICK, color="white", zorder=5)
                else:
                    target_ax.text(val + 0.3, idx, f"{int(reads):,} reads", va="center", ha="left", fontsize=FONT_YTICK, color="black", zorder=5)

        parent_label_name = parent_col.replace("_name", "").replace("_", " ").title()
        handles = [mpatches.Patch(facecolor=color_map[g], edgecolor="none", label=g) for g in unique_groups]
        ax.legend(handles=handles, title=parent_label_name, fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE, bbox_to_anchor=(1.01 if not axis_break else 1.25, 0.5), loc="center left", frameon=False)

        ax.set_xlabel(f"% positive samples (n = {n_tot})", fontsize=FONT_XLABEL)
        level_title = "Accession Number" if level_label == "Accession" else level_label
        ax.set_title(f"Distribution of Viruses – {level_title} Level", fontsize=FONT_TITLE)

        if axis_break:
            ax.set_xlim(0, x_cut_low)
            ax_right.set_xlim(x_cut_high, x_max_right)
            for a in axes_list:
                a.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
            draw_break_markers(ax, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
        else:
            ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
            ax.set_xlim(0, max_val * 1.35)

        for a in axes_list:
            a.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
            a.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

        fig.text(0.5, 0.005, f"source: {strat_tag}", ha="center", fontsize=FONT_FOOTNOTE, fontstyle="normal", color="black")
        plt.tight_layout(rect=[0, 0.04, 1, 1])
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig)
        print(f"  ✓ saved {save_path}")


def plot_accession_family_overlay(data, save_path_prefix, n_total_cohort=0, strat_tag="", enable_broken=True):
    """
    Accession plot with family-level translucent overlay background boxes.
    Generates both brokenaxis and nonbrokenaxis versions.
    """
    above = ["family_name", "genus_name", "species_name"]
    group_cols = ["virus_accession", "virus_name_sanitized"] + above

    grp = (
        data.groupby(group_cols)
        .agg(
            n_positive_samples=("sample_id", "nunique"),
            total_reads=("virus_mapped_reads", "sum"),
        )
        .reset_index()
    )
    if grp.empty:
        return

    n_tot = n_total_cohort if n_total_cohort > 0 else data["sample_id"].nunique()
    grp["pct_positive"] = grp["n_positive_samples"] / float(n_tot) * 100.0 if n_tot > 0 else 0.0
    grp["label"] = grp["virus_accession"].astype(str) + ": " + grp["virus_name_sanitized"].astype(str)

    sort_cols = []
    asc_flags = []
    for c in above:
        c_tot = data.groupby(c)["sample_id"].nunique().to_dict()
        grp[c + "_tot"] = grp[c].map(c_tot)
        grp[c + "_unk"] = grp[c] == "Unknown"
        sort_cols.extend([c + "_unk", c + "_tot", c])
        asc_flags.extend([True, False, True])

    sort_cols.extend(["n_positive_samples", "total_reads"])
    asc_flags.extend([False, False])
    grp = grp.sort_values(by=sort_cols, ascending=asc_flags).reset_index(drop=True)

    genus_order = list(dict.fromkeys(grp["genus_name"]))
    genus_color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(genus_order)}
    acc_bar_colors = [genus_color_map[g] for g in grp["genus_name"]]

    fam_order = list(dict.fromkeys(grp["family_name"]))
    fam_samples = data.groupby("family_name")["sample_id"].nunique().to_dict()
    fam_reads = data.groupby("family_name")["virus_mapped_reads"].sum().to_dict()

    fam_palette = ["#C6D8EF", "#FEE8C8", "#E7D4E8", "#FDDBC7", "#E0ECF4", "#E5F5E0"]
    fam_tint_map = {f: fam_palette[i % len(fam_palette)] for i, f in enumerate(fam_order)}

    fam_pcts = [v / float(n_tot) * 100.0 for v in fam_samples.values()]
    all_vals = list(grp["pct_positive"]) + fam_pcts
    axis_break_detected = detect_axis_break(all_vals, "samples") if enable_broken else None

    break_variants = []
    if axis_break_detected:
        break_variants.append((axis_break_detected, f"{save_path_prefix}__brokenaxis.tiff"))
    break_variants.append((None, f"{save_path_prefix}__nonbrokenaxis.tiff"))

    n_bars = len(grp)
    fig_h = max(4.8, n_bars * 0.40)
    y_positions = list(range(len(grp)))

    for axis_break, save_path in break_variants:
        if axis_break:
            x_cut_low, x_cut_high, x_max_right = axis_break
            fig, (ax, ax_right) = plt.subplots(1, 2, figsize=(12.5, fig_h), sharey=True, gridspec_kw={"width_ratios": [3.2, 1.0], "wspace": 0.045})
            axes_list = [ax, ax_right]
            for a in axes_list:
                a.spines["top"].set_visible(False)
                a.spines["right"].set_visible(False)
                a.spines["left"].set_linewidth(0.6)
                a.spines["bottom"].set_linewidth(0.6)
            ax_right.spines["left"].set_visible(False)
            ax_right.tick_params(left=False)
        else:
            fig, ax = plt.subplots(figsize=(11.5, fig_h))
            axes_list = [ax]
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(0.6)
            ax.spines["bottom"].set_linewidth(0.6)

        crossing_bars = []
        for fam in fam_order:
            idx_list = grp.index[grp["family_name"] == fam].tolist()
            if not idx_list:
                continue
            y_min = min(idx_list) - FAM_BOX_MARGIN
            y_max = max(idx_list) + FAM_BOX_MARGIN
            height = y_max - y_min
            pct = fam_samples[fam] / float(n_tot) * 100.0
            tint = fam_tint_map[fam]

            for a in axes_list:
                rect = plt.Rectangle((0, y_min), pct, height, facecolor=tint, alpha=0.55, edgecolor="none", zorder=1)
                a.add_patch(rect)

            fam_short = fam.removesuffix("viridae")
            fam_label = f"{fam_short}: n={fam_samples[fam]} ({pct:.1f}%), {int(fam_reads[fam]):,} reads"
            if axis_break and pct > x_cut_low:
                ax_right.text(pct + 0.5, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=FONT_FAMILY, color="#1A252F", zorder=2)
            else:
                x_pos = pct + 1.0 if pct > 10 else 13.0
                ax.text(x_pos, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=FONT_FAMILY, color="#1A252F", zorder=2)

        for a in axes_list:
            bars = a.barh(y_positions, grp["pct_positive"], height=BAR_HEIGHT, color=acc_bar_colors, edgecolor="none", zorder=3)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(grp["label"], fontsize=FONT_YTICK)
        ax.invert_yaxis()

        for idx, (w, reads) in enumerate(zip(grp["pct_positive"], grp["total_reads"])):
            if axis_break and w > x_cut_low:
                crossing_bars.append(idx)
                ax_right.text(w + 0.3, idx, f"{int(reads):,} reads", va="center", ha="left", fontsize=FONT_YTICK, color="black", zorder=4)
            else:
                if w > 12 and not axis_break:
                    ax.text(w - 1, idx, f"{int(reads):,} reads", va="center", ha="right", fontsize=FONT_YTICK, color="white", zorder=4)
                else:
                    ax.text(w + 0.25, idx, f"{int(reads):,} reads", va="center", ha="left", fontsize=FONT_YTICK, color="black", zorder=4)

        handles = [mpatches.Patch(facecolor=genus_color_map[g], edgecolor="none", label=g) for g in genus_order]
        ax.legend(handles=handles, title="Genus", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE, bbox_to_anchor=(1.01 if not axis_break else 1.25, 0.5), loc="center left", frameon=False)

        ax.set_xlabel(f"% positive samples (n = {n_tot})", fontsize=FONT_XLABEL)
        ax.set_title("Distribution of Viruses – Accession Number Level with Family Overlay", fontsize=FONT_TITLE)

        if axis_break:
            ax.set_xlim(0, x_cut_low)
            ax_right.set_xlim(x_cut_high, x_max_right)
            for a in axes_list:
                a.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
            draw_break_markers(ax, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
        else:
            max_fam_val = max(fam_samples.values()) / float(n_tot) * 100.0 if fam_samples else 1.0
            ax.set_xlim(0, max_fam_val * 1.55)
            ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))

        for a in axes_list:
            a.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
            a.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

        fig.text(0.5, 0.005, f"source: {strat_tag}   •   family overlay names shown without '-viridae' suffix", ha="center", fontsize=FONT_FOOTNOTE, fontstyle="normal", color="black")
        plt.tight_layout(rect=[0, 0.04, 1, 1])
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig)
        print(f"  ✓ saved {save_path}")


def render_sample_enrichment_overlay(
    df_sub,
    out_path,
    title_str,
    n_total_cohort=0,
    strat_tag=SOURCE_FILE_FILTER,
    mode="samples",
):
    """
    Renders a 300 DPI publication horizontal sample-enrichment bar plot:
    - Grouped by viral family with family-level shaded background boxes.
    - Horizontal bars stacked/segmented by sample ID in distinct colors with chunky bar height (0.88).
    - Top and right spines removed.
    - Standardized typography and natural canvas height scaling across all cohorts.
    """
    if df_sub.empty or "family_name" not in df_sub.columns or "genus_name" not in df_sub.columns:
        return

    above_cols = ["family_name", "genus_name", "species_name"]
    grp_cols = ["virus_accession", "virus_name_sanitized"] + [c for c in above_cols if c in df_sub.columns]

    grp_data = df_sub.groupby(grp_cols).agg(
        n_positive_samples=("sample_id", "nunique"),
        total_reads=("virus_mapped_reads", "sum"),
        total_copy_number=("viral_copy_number", "sum"),
    ).reset_index()

    if mode == "samples":
        grp_data["pct_positive"] = grp_data["n_positive_samples"] / float(n_total_cohort) * 100.0 if n_total_cohort > 0 else 0.0
    grp_data["label"] = grp_data["virus_accession"].astype(str) + ": " + grp_data["virus_name_sanitized"].astype(str)

    # Hierarchical sorting
    s_cols_sub, s_asc_sub = [], []
    for c in above_cols:
        if c in df_sub.columns:
            if mode == "samples":
                c_tot = df_sub.groupby(c)["sample_id"].nunique().to_dict()
            elif mode == "copynumber":
                c_tot = df_sub.groupby(c)["viral_copy_number"].sum().to_dict()
            else:
                c_tot = df_sub.groupby(c)["virus_mapped_reads"].sum().to_dict()
            grp_data[c + "_tot"] = grp_data[c].map(c_tot)
            grp_data[c + "_unk"] = grp_data[c] == "Unknown"
            s_cols_sub.extend([c + "_unk", c + "_tot", c])
            s_asc_sub.extend([True, False, True])

    if mode == "samples":
        s_cols_sub.extend(["n_positive_samples", "total_reads"])
    elif mode == "copynumber":
        s_cols_sub.extend(["total_copy_number", "total_reads", "n_positive_samples"])
    else:
        s_cols_sub.extend(["total_reads", "n_positive_samples"])
    s_asc_sub.extend([False] * (len(s_cols_sub) - len(s_asc_sub)))
    grp_data = grp_data.sort_values(by=s_cols_sub, ascending=s_asc_sub).reset_index(drop=True)

    if grp_data.empty:
        return

    # Assign a distinct color to each sample
    unique_samples_sub = list(dict.fromkeys(df_sub["sample_id"]))
    sample_color_map = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(unique_samples_sub)}

    # Family background tinting
    fam_order_sub = list(dict.fromkeys(grp_data["family_name"]))
    fam_samples_sub = df_sub.groupby("family_name")["sample_id"].nunique().to_dict()
    fam_reads_sub = df_sub.groupby("family_name")["virus_mapped_reads"].sum().to_dict()
    fam_cn_sub = df_sub.groupby("family_name")["viral_copy_number"].sum().to_dict()

    fam_palette = ["#C6D8EF", "#FEE8C8", "#E7D4E8", "#FDDBC7", "#E0ECF4", "#E5F5E0"]
    fam_tint_map = {f: fam_palette[i % len(fam_palette)] for i, f in enumerate(fam_order_sub)}

    # Determine if an outlier gap exists for broken axis
    if mode == "samples":
        pct_fam = [v / float(n_total_cohort) * 100.0 for v in fam_samples_sub.values()]
        all_vals = list(grp_data["pct_positive"]) + pct_fam if n_total_cohort > 0 else []
    elif mode == "copynumber":
        all_vals = list(grp_data["total_copy_number"]) + list(fam_cn_sub.values())
    else:
        all_vals = list(grp_data["total_reads"]) + list(fam_reads_sub.values())

    axis_break_detected = detect_axis_break(all_vals, mode)

    n_unique_samples = len(unique_samples_sub)
    if n_unique_samples > 100:
        ncol_legend = 10
    elif n_unique_samples > 40:
        ncol_legend = 8
    else:
        ncol_legend = min(6, max(1, n_unique_samples))

    n_rows = int(np.ceil(n_unique_samples / ncol_legend))
    plot_h = max(4.8, len(grp_data) * 0.40)
    legend_h = n_rows * 0.18 + 0.6
    total_h = plot_h + legend_h

    y_pos_sub = list(range(len(grp_data)))

    # Compute max value for x-axis scaling
    if mode == "samples":
        max_fam_val = max(fam_samples_sub.values()) / float(n_total_cohort) * 100.0 if fam_samples_sub and n_total_cohort > 0 else 10.0
    elif mode == "copynumber":
        max_fam_val = max(fam_cn_sub.values()) if fam_cn_sub else 0.01
    else:
        max_fam_val = max(fam_reads_sub.values()) if fam_reads_sub else 10.0

    # Determine output path variants:
    # If break detected: generate both '__brokenaxis' and '__nonbrokenaxis' versions.
    # Otherwise: generate '__nonbrokenaxis' version.
    parent = out_path.parent
    stem = out_path.stem
    ext = out_path.suffix
    if axis_break_detected:
        break_variants = [
            (axis_break_detected, parent / f"{stem}__brokenaxis{ext}"),
            (None,                parent / f"{stem}__nonbrokenaxis{ext}"),
        ]
    else:
        break_variants = [
            (None,                parent / f"{stem}__nonbrokenaxis{ext}"),
        ]

    for axis_break, save_path in break_variants:
        if axis_break:
            x_cut_low, x_cut_high, x_max_right = axis_break
            fig_sub, (ax_sub, ax_right) = plt.subplots(
                1, 2, sharey=True, figsize=(11.5, total_h), dpi=300,
                gridspec_kw={"width_ratios": [3.0, 1.1], "wspace": 0.04}
            )
            axes_list = [ax_sub, ax_right]
        else:
            fig_sub, ax_sub = plt.subplots(figsize=(11.5, total_h), dpi=300)
            ax_right = None
            axes_list = [ax_sub]

        bottom_frac = (legend_h + 0.3) / total_h
        top_frac = 1.0 - (0.45 / total_h)
        fig_sub.subplots_adjust(bottom=bottom_frac, top=top_frac, left=0.22, right=0.95)

        # Spine styling
        for a in axes_list:
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
            a.spines["bottom"].set_linewidth(0.6)
        ax_sub.spines["left"].set_linewidth(0.6)
        if ax_right:
            ax_right.spines["left"].set_visible(False)
            ax_right.yaxis.set_ticks_position("none")
            ax_right.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
            ax_right.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)

        # 1. Draw Family shaded background boxes (framing rows with FAM_BOX_MARGIN)
        for fam in fam_order_sub:
            idx_list = grp_data.index[grp_data["family_name"] == fam].tolist()
            if not idx_list:
                continue
            y_min = min(idx_list) - FAM_BOX_MARGIN
            y_max = max(idx_list) + FAM_BOX_MARGIN
            tint = fam_tint_map[fam]
            fam_short = fam.removesuffix("viridae")

            if mode == "samples":
                val = fam_samples_sub[fam] / float(n_total_cohort) * 100.0 if n_total_cohort > 0 else 0.0
                fam_label = f"{fam_short}: n={fam_samples_sub[fam]} ({val:.1f}%), {int(fam_reads_sub[fam]):,} reads"
            elif mode == "copynumber":
                val = fam_cn_sub[fam]
                fam_label = f"{fam_short}: {val:.4g} copies (n={fam_samples_sub[fam]})"
            else:
                val = fam_reads_sub[fam]
                fam_label = f"{fam_short}: {int(val):,} reads (n={fam_samples_sub[fam]})"

            for a in axes_list:
                rect = plt.Rectangle((0, y_min), val, y_max - y_min, facecolor=tint, alpha=0.55, edgecolor="none", zorder=1)
                a.add_patch(rect)

            if axis_break and val >= x_cut_high:
                x_pos = val + (x_max_right - x_cut_high) * 0.02
                ax_right.text(x_pos, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=FONT_FAMILY, color="#1A252F", zorder=2)
            else:
                ref_max = x_cut_low if axis_break else max_fam_val
                if mode == "copynumber":
                    x_pos = val + ref_max * 0.02
                else:
                    x_pos = val + max(ref_max * 0.02, 1.0)
                ax_sub.text(x_pos, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=FONT_FAMILY, color="#1A252F", zorder=2)

        # 2. Draw segmented horizontal bars stacked per sample (chunky height=0.88)
        crossing_bars = []
        for idx_row, row_v in grp_data.iterrows():
            v_acc = row_v["virus_accession"]
            v_df = df_sub[df_sub["virus_accession"] == v_acc].drop_duplicates(subset=["sample_id"])
            left_val = 0.0
            if axis_break:
                if mode == "samples":
                    row_tot = float(row_v.get("pct_positive", 0.0))
                elif mode == "copynumber":
                    row_tot = float(row_v.get("total_copy_number", 0.0))
                else:
                    row_tot = float(row_v.get("total_reads", 0.0))
                if row_tot >= x_cut_high:
                    crossing_bars.append(idx_row)
            for s_row in v_df.itertuples(index=False):
                s_id = s_row.sample_id
                if mode == "samples":
                    seg_w = 1.0 / float(n_total_cohort) * 100.0 if n_total_cohort > 0 else 0.0
                elif mode == "copynumber":
                    seg_w = float(getattr(s_row, "viral_copy_number", 0.0))
                else:
                    seg_w = float(s_row.virus_mapped_reads)
                if seg_w > 0:
                    for a in axes_list:
                        a.barh(idx_row, seg_w, left=left_val, height=BAR_HEIGHT, color=sample_color_map[s_id], edgecolor="none", zorder=3)
                    left_val += seg_w

        ax_sub.set_yticks(y_pos_sub)
        ax_sub.set_yticklabels(grp_data["label"], fontsize=FONT_YTICK)
        ax_sub.invert_yaxis()

        # 3. Position sample color legend UNDER the plot below the X-axis
        sample_handles = [mpatches.Patch(facecolor=sample_color_map[s], edgecolor="none", label=s) for s in unique_samples_sub]
        ax_sub.legend(
            handles=sample_handles, title="Samples", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE,
            bbox_to_anchor=(0.64 if axis_break else 0.50, -0.10), loc="upper center", ncol=ncol_legend, frameon=False,
        )

        if mode == "samples":
            ax_sub.set_xlabel(f"% positive samples (n = {n_total_cohort})", fontsize=FONT_XLABEL, x=0.64 if axis_break else 0.50)
            if axis_break:
                ax_sub.set_xlim(0, x_cut_low)
                ax_right.set_xlim(x_cut_high, x_max_right)
                for a in axes_list:
                    a.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
                draw_break_markers(ax_sub, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
            else:
                ax_sub.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
                ax_sub.set_xlim(0, max_fam_val * 1.55)
        elif mode == "copynumber":
            ax_sub.set_xlabel("Viral copy number", fontsize=FONT_XLABEL, x=0.64 if axis_break else 0.50)
            if axis_break:
                ax_sub.set_xlim(0, x_cut_low)
                ax_right.set_xlim(x_cut_high, x_max_right)
                for a in axes_list:
                    cn_fmt = mtick.ScalarFormatter(useOffset=False)
                    cn_fmt.set_scientific(False)
                    a.xaxis.set_major_formatter(cn_fmt)
                draw_break_markers(ax_sub, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
            else:
                cn_fmt = mtick.ScalarFormatter(useOffset=False)
                cn_fmt.set_scientific(False)
                ax_sub.xaxis.set_major_formatter(cn_fmt)
                ax_sub.set_xlim(0, max_fam_val * 1.55)
        else:
            ax_sub.set_xlabel("Number of mapped reads", fontsize=FONT_XLABEL, x=0.64 if axis_break else 0.50)
            if axis_break:
                ax_sub.set_xlim(0, x_cut_low)
                ax_right.set_xlim(x_cut_high, x_max_right)
                for a in axes_list:
                    a.xaxis.set_major_formatter(mtick.StrMethodFormatter("{x:,.0f}"))
                draw_break_markers(ax_sub, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
            else:
                ax_sub.xaxis.set_major_formatter(mtick.StrMethodFormatter("{x:,.0f}"))
                ax_sub.set_xlim(0, max(max_fam_val * 1.45, max_fam_val + 10.0))

        ax_sub.set_title(title_str, fontsize=FONT_TITLE, pad=8, x=0.64 if axis_break else 0.50)
        ax_sub.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
        ax_sub.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

        fig_sub.text(0.5, 0.01, f"source: {strat_tag}   •   family overlay names shown without '-viridae' suffix", ha="center", fontsize=FONT_FOOTNOTE, color="#333333")

        fig_sub.savefig(save_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig_sub)
        print(f"  ✓ saved {save_path}")


# ── Helper: specimen class overlay plot (full cohort) ────────────────────────
CLASS_COLORS = {
    "Interneurons (INT)": "#0072B2",
    "Astrocytes & Microglia (ASTMIG)": "#E69F00",
    "Oligodendrocytes (OLI)": "#009E73",
    "Pyramidal neurons (PYR)": "#D55E00",
    "Medium spiny neurons (MSN)": "#7570B3",
    "Bulk tissue": "#666666",
}
CLASS_ORDER = list(CLASS_COLORS.keys())

def get_specimen_class(sample_id):
    s = str(sample_id)
    if "INT" in s: return "Interneurons (INT)"
    if "ASTMIG" in s: return "Astrocytes & Microglia (ASTMIG)"
    if "OLI" in s: return "Oligodendrocytes (OLI)"
    if "PYR" in s: return "Pyramidal neurons (PYR)"
    if "MSN" in s: return "Medium spiny neurons (MSN)"
    return "Bulk tissue"

def render_specimen_class_overlay(
    df_cohort,
    out_path,
    title_str,
    n_total_cohort=N_TOTAL,
    strat_tag=SOURCE_FILE_FILTER,
    mode="samples",
):
    """
    Renders a 300 DPI publication horizontal sample-enrichment bar plot for the full cohort,
    segmenting each virus bar by 6 biologically distinct specimen classes rather than 241 individual sample IDs.
    """
    if df_cohort.empty:
        return

    df_plot = df_cohort.copy()
    df_plot["specimen_class"] = df_plot["sample_id"].apply(get_specimen_class)

    above_cols = ["family_name", "genus_name", "species_name"]
    grp_cols = ["virus_accession", "virus_name_sanitized"] + [c for c in above_cols if c in df_plot.columns]

    grp_data = df_plot.groupby(grp_cols).agg(
        n_positive_samples=("sample_id", "nunique"),
        total_reads=("virus_mapped_reads", "sum"),
        total_copy_number=("viral_copy_number", "sum"),
    ).reset_index()

    if mode == "samples":
        grp_data["pct_positive"] = grp_data["n_positive_samples"] / float(n_total_cohort) * 100.0 if n_total_cohort > 0 else 0.0
    grp_data["label"] = grp_data["virus_accession"].astype(str) + ": " + grp_data["virus_name_sanitized"].astype(str)

    # Hierarchical sorting
    s_cols_sub, s_asc_sub = [], []
    for c in above_cols:
        if c in df_plot.columns:
            if mode == "samples":
                c_tot = df_plot.groupby(c)["sample_id"].nunique().to_dict()
            elif mode == "copynumber":
                c_tot = df_plot.groupby(c)["viral_copy_number"].sum().to_dict()
            else:
                c_tot = df_plot.groupby(c)["virus_mapped_reads"].sum().to_dict()
            grp_data[c + "_tot"] = grp_data[c].map(c_tot)
            grp_data[c + "_unk"] = grp_data[c] == "Unknown"
            s_cols_sub.extend([c + "_unk", c + "_tot", c])
            s_asc_sub.extend([True, False, True])

    if mode == "samples":
        s_cols_sub.extend(["n_positive_samples", "total_reads"])
    elif mode == "copynumber":
        s_cols_sub.extend(["total_copy_number", "total_reads", "n_positive_samples"])
    else:
        s_cols_sub.extend(["total_reads", "n_positive_samples"])
    s_asc_sub.extend([False] * (len(s_cols_sub) - len(s_asc_sub)))
    grp_data = grp_data.sort_values(by=s_cols_sub, ascending=s_asc_sub).reset_index(drop=True)

    fam_order_sub = list(dict.fromkeys(grp_data["family_name"]))
    fam_samples_sub = df_plot.groupby("family_name")["sample_id"].nunique().to_dict()
    fam_reads_sub = df_plot.groupby("family_name")["virus_mapped_reads"].sum().to_dict()
    fam_cn_sub = df_plot.groupby("family_name")["viral_copy_number"].sum().to_dict()

    fam_palette = ["#C6D8EF", "#FEE8C8", "#E7D4E8", "#FDDBC7", "#E0ECF4", "#E5F5E0"]
    fam_tint_map = {f: fam_palette[i % len(fam_palette)] for i, f in enumerate(fam_order_sub)}

    # Determine if an outlier gap exists for broken axis
    if mode == "samples":
        pct_fam = [v / float(n_total_cohort) * 100.0 for v in fam_samples_sub.values()]
        all_vals = list(grp_data["pct_positive"]) + pct_fam if n_total_cohort > 0 else []
    elif mode == "copynumber":
        all_vals = list(grp_data["total_copy_number"]) + list(fam_cn_sub.values())
    else:
        all_vals = list(grp_data["total_reads"]) + list(fam_reads_sub.values())

    axis_break_detected = detect_axis_break(all_vals, mode)

    plot_h = max(4.8, len(grp_data) * 0.40)
    legend_h = 1.0
    total_h = plot_h + legend_h

    y_pos = list(range(len(grp_data)))

    if mode == "samples":
        max_fam_val = max(fam_samples_sub.values()) / float(n_total_cohort) * 100.0 if fam_samples_sub and n_total_cohort > 0 else 10.0
    elif mode == "copynumber":
        max_fam_val = max(fam_cn_sub.values()) if fam_cn_sub else 0.01
    else:
        max_fam_val = max(fam_reads_sub.values()) if fam_reads_sub else 10.0

    # Determine output path variants:
    # If break detected: generate both '__brokenaxis' and '__nonbrokenaxis' versions.
    # Otherwise: generate '__nonbrokenaxis' version.
    parent = out_path.parent
    stem = out_path.stem
    ext = out_path.suffix
    if axis_break_detected:
        break_variants = [
            (axis_break_detected, parent / f"{stem}__brokenaxis{ext}"),
            (None,                parent / f"{stem}__nonbrokenaxis{ext}"),
        ]
    else:
        break_variants = [
            (None,                parent / f"{stem}__nonbrokenaxis{ext}"),
        ]

    for axis_break, save_path in break_variants:
        if axis_break:
            x_cut_low, x_cut_high, x_max_right = axis_break
            fig, (ax, ax_right) = plt.subplots(
                1, 2, sharey=True, figsize=(11.5, total_h), dpi=300,
                gridspec_kw={"width_ratios": [3.0, 1.1], "wspace": 0.04}
            )
            axes_list = [ax, ax_right]
        else:
            fig, ax = plt.subplots(figsize=(11.5, total_h), dpi=300)
            ax_right = None
            axes_list = [ax]

        bottom_frac = (legend_h + 0.3) / total_h
        top_frac = 1.0 - (0.45 / total_h)
        fig.subplots_adjust(bottom=bottom_frac, top=top_frac, left=0.22, right=0.95)

        # Spine styling
        for a in axes_list:
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
            a.spines["bottom"].set_linewidth(0.6)
        ax.spines["left"].set_linewidth(0.6)
        if ax_right:
            ax_right.spines["left"].set_visible(False)
            ax_right.yaxis.set_ticks_position("none")
            ax_right.tick_params(axis="y", which="both", left=False, right=False, labelleft=False, labelright=False)
            ax_right.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)

        # Family boxes
        for fam in fam_order_sub:
            idx_list = grp_data.index[grp_data["family_name"] == fam].tolist()
            if not idx_list:
                continue
            y_min = min(idx_list) - FAM_BOX_MARGIN
            y_max = max(idx_list) + FAM_BOX_MARGIN
            tint = fam_tint_map[fam]
            fam_short = fam.removesuffix("viridae")

            if mode == "samples":
                val = fam_samples_sub[fam] / float(n_total_cohort) * 100.0 if n_total_cohort > 0 else 0.0
                fam_label = f"{fam_short}: n={fam_samples_sub[fam]} ({val:.1f}%), {int(fam_reads_sub[fam]):,} reads"
            elif mode == "copynumber":
                val = fam_cn_sub[fam]
                fam_label = f"{fam_short}: {val:.4g} copies (n={fam_samples_sub[fam]})"
            else:
                val = fam_reads_sub[fam]
                fam_label = f"{fam_short}: {int(val):,} reads (n={fam_samples_sub[fam]})"

            for a in axes_list:
                rect = plt.Rectangle((0, y_min), val, y_max - y_min, facecolor=tint, alpha=0.55, edgecolor="none", zorder=1)
                a.add_patch(rect)

            if axis_break and val >= x_cut_high:
                x_pos = val + (x_max_right - x_cut_high) * 0.02
                ax_right.text(x_pos, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=FONT_FAMILY, color="#1A252F", zorder=2)
            else:
                ref_max = x_cut_low if axis_break else max_fam_val
                if mode == "copynumber":
                    x_pos = val + ref_max * 0.02
                else:
                    x_pos = val + max(ref_max * 0.02, 1.0)
                ax.text(x_pos, (y_min + y_max) / 2, fam_label, va="center", ha="left", fontsize=FONT_FAMILY, color="#1A252F", zorder=2)

        # Bars by specimen class
        crossing_bars = []
        for idx_row, row_v in grp_data.iterrows():
            v_acc = row_v["virus_accession"]
            v_df = df_plot[df_plot["virus_accession"] == v_acc]
            left_val = 0.0
            if axis_break:
                if mode == "samples":
                    row_tot = float(row_v.get("pct_positive", 0.0))
                elif mode == "copynumber":
                    row_tot = float(row_v.get("total_copy_number", 0.0))
                else:
                    row_tot = float(row_v.get("total_reads", 0.0))
                if row_tot >= x_cut_high:
                    crossing_bars.append(idx_row)
            for cls_name in CLASS_ORDER:
                cls_sub = v_df[v_df["specimen_class"] == cls_name]
                if cls_sub.empty:
                    continue
                if mode == "samples":
                    seg_w = cls_sub["sample_id"].nunique() / float(n_total_cohort) * 100.0
                elif mode == "copynumber":
                    seg_w = float(cls_sub["viral_copy_number"].sum())
                else:
                    seg_w = float(cls_sub["virus_mapped_reads"].sum())
                if seg_w > 0:
                    for a in axes_list:
                        a.barh(idx_row, seg_w, left=left_val, height=BAR_HEIGHT, color=CLASS_COLORS[cls_name], edgecolor="none", zorder=3)
                    left_val += seg_w

        ax.set_yticks(y_pos)
        ax.set_yticklabels(grp_data["label"], fontsize=FONT_YTICK)
        ax.invert_yaxis()

        handles = [mpatches.Patch(facecolor=CLASS_COLORS[c], edgecolor="none", label=c) for c in CLASS_ORDER]
        ax.legend(
            handles=handles, title="Specimen class", fontsize=FONT_LEGEND, title_fontsize=FONT_LEGEND_TITLE,
            bbox_to_anchor=(0.64 if axis_break else 0.50, -0.09), loc="upper center", ncol=3, frameon=False, columnspacing=1.0,
        )

        if mode == "samples":
            ax.set_xlabel(f"% positive samples (n = {n_total_cohort})", fontsize=FONT_XLABEL, x=0.64 if axis_break else 0.50)
            if axis_break:
                ax.set_xlim(0, x_cut_low)
                ax_right.set_xlim(x_cut_high, x_max_right)
                for a in axes_list:
                    a.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
                draw_break_markers(ax, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
            else:
                ax.xaxis.set_major_formatter(mtick.PercentFormatter(decimals=1))
                ax.set_xlim(0, max_fam_val * 1.55)
        elif mode == "copynumber":
            ax.set_xlabel("Viral copy number", fontsize=FONT_XLABEL, x=0.64 if axis_break else 0.50)
            if axis_break:
                ax.set_xlim(0, x_cut_low)
                ax_right.set_xlim(x_cut_high, x_max_right)
                for a in axes_list:
                    cn_fmt = mtick.ScalarFormatter(useOffset=False)
                    cn_fmt.set_scientific(False)
                    a.xaxis.set_major_formatter(cn_fmt)
                draw_break_markers(ax, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
            else:
                cn_fmt = mtick.ScalarFormatter(useOffset=False)
                cn_fmt.set_scientific(False)
                ax.xaxis.set_major_formatter(cn_fmt)
                ax.set_xlim(0, max_fam_val * 1.55)
        else:
            ax.set_xlabel("Number of mapped reads", fontsize=FONT_XLABEL, x=0.64 if axis_break else 0.50)
            if axis_break:
                ax.set_xlim(0, x_cut_low)
                ax_right.set_xlim(x_cut_high, x_max_right)
                for a in axes_list:
                    a.xaxis.set_major_formatter(mtick.StrMethodFormatter("{x:,.0f}"))
                draw_break_markers(ax, ax_right, crossing_bars=crossing_bars, bar_height=BAR_HEIGHT)
            else:
                ax.xaxis.set_major_formatter(mtick.StrMethodFormatter("{x:,.0f}"))
                ax.set_xlim(0, max(max_fam_val * 1.45, max_fam_val + 10.0))

        ax.set_title(title_str, fontsize=FONT_TITLE, pad=8, x=0.64 if axis_break else 0.50)
        ax.tick_params(axis="x", labelsize=FONT_XTICK, direction="out", length=3, width=0.6)
        ax.tick_params(axis="y", labelsize=FONT_YTICK, direction="out", length=3, width=0.6)

        fig.text(0.5, 0.01, f"source: {strat_tag}   •   family overlay names shown without '-viridae' suffix", ha="center", fontsize=FONT_FOOTNOTE, color="#333333")
        fig.savefig(save_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
        plt.close(fig)
        print(f"  ✓ saved {save_path}")


# ── Generate taxonomy plots ──────────────────────────────────────────────────
taxonomy_levels = [
    ("virus_accession",  "virus_name_sanitized", "Accession",  "taxonomy_accession.tiff"),
    ("species_taxid",    "species_name",         "Species",    "taxonomy_species.tiff"),
    ("genus_taxid",      "genus_name",           "Genus",      "taxonomy_genus.tiff"),
    ("family_taxid",     "family_name",           "Family",     "taxonomy_family.tiff"),
    ("order_taxid",      "order_name",            "Order",      "taxonomy_order.tiff"),
    ("realm_taxid",      "realm_name",            "Realm",      "taxonomy_realm.tiff"),
]

print("\n── Taxonomy plots ──")
for taxid_col, name_col, level_label, fname in taxonomy_levels:
    plot_taxonomy_level(df_filt, taxid_col, name_col, level_label, fname)

# Hybrid accession plot with family overlay
plot_accession_family_overlay(df_filt, "taxonomy_accession_family_overlay.tiff")

# ── Sample enrichment & read abundance overlay plots ─────────────────────────
overlay_modes = [
    ("samples", "Sample Enrichment", "sample_enrichment_overlay_samples"),
    ("reads",   "Read Abundance",    "sample_enrichment_overlay_reads"),
    ("copynumber", "Viral Copy Number", "sample_enrichment_overlay_copynumber"),
]

if BSMN_LABELS_PATH.exists():
    with open(BSMN_LABELS_PATH, "r") as f:
        bsmn_labels = [line.strip() for line in f if line.strip()]

    # Pre-filter labeled and unlabeled datasets
    sub_labeled = df_filt[df_filt["sample_id"].apply(lambda s: any(l in str(s) for l in bsmn_labels))].copy()
    n_tot_labeled = df[df["sample_id"].apply(lambda s: any(l in str(s) for l in bsmn_labels))]["sample_id"].nunique()

    sub_unlabeled = df_filt[~df_filt["sample_id"].apply(lambda s: any(l in str(s) for l in bsmn_labels))].copy()
    n_tot_unlabeled = df[~df["sample_id"].apply(lambda s: any(l in str(s) for l in bsmn_labels))]["sample_id"].nunique()

    for mode, mode_title, file_prefix in overlay_modes:
        print(f"\n── Overlay plots: {mode_title} (BSMN sorted labels & cohorts) ──")

        # Individual BSMN sorted labels (exact 1800 x 1800 px square collage panels)
        for lbl in bsmn_labels:
            sub_df = df_filt[df_filt["sample_id"].str.contains(lbl, regex=False)].copy()
            if sub_df.empty:
                continue

            n_label_cohort = df[df["sample_id"].str.contains(lbl, regex=False)]["sample_id"].nunique()
            if n_label_cohort == 0:
                n_label_cohort = sub_df["sample_id"].nunique()

            out_fname = f"{file_prefix}_{lbl}.tiff"
            render_sample_enrichment_overlay(
                df_sub=sub_df,
                out_path=OUT_DIR / out_fname,
                title_str=f"Distribution of Viruses – {mode_title} ({lbl})",
                n_total_cohort=n_label_cohort,
                strat_tag=SOURCE_FILE_FILTER,
                mode=mode,
            )

        # 1. All labeled samples
        if not sub_labeled.empty:
            out_fname = f"{file_prefix}_all_labeled.tiff"
            render_sample_enrichment_overlay(
                df_sub=sub_labeled,
                out_path=OUT_DIR / out_fname,
                title_str=f"Distribution of Viruses – {mode_title} Overlay (All Labeled Samples)",
                n_total_cohort=n_tot_labeled,
                strat_tag=SOURCE_FILE_FILTER,
                mode=mode,
            )

        # 2. Unlabeled samples
        if not sub_unlabeled.empty:
            out_fname = f"{file_prefix}_unlabeled.tiff"
            render_sample_enrichment_overlay(
                df_sub=sub_unlabeled,
                out_path=OUT_DIR / out_fname,
                title_str=f"Distribution of Viruses – {mode_title} Overlay (Unlabeled Samples)",
                n_total_cohort=n_tot_unlabeled,
                strat_tag=SOURCE_FILE_FILTER,
                mode=mode,
            )

        # 3. All samples (full cohort) - segmented by individual sample ID
        out_fname = f"{file_prefix}_all_samples.tiff"
        render_sample_enrichment_overlay(
            df_sub=df_filt.copy(),
            out_path=OUT_DIR / out_fname,
            title_str=f"Distribution of Viruses – {mode_title} Overlay (All Samples)",
            n_total_cohort=N_TOTAL,
            strat_tag=SOURCE_FILE_FILTER,
            mode=mode,
        )

        # 4. All samples (full cohort) - segmented by specimen class (NEW ADDITIONAL PLOT)
        out_fname_class = f"{file_prefix}_all_samples_by_class.tiff"
        render_specimen_class_overlay(
            df_cohort=df_filt.copy(),
            out_path=OUT_DIR / out_fname_class,
            title_str=f"Distribution of Viruses – {mode_title} by Specimen Class (All Samples)",
            n_total_cohort=N_TOTAL,
            strat_tag=SOURCE_FILE_FILTER,
            mode=mode,
        )
else:
    for mode, mode_title, file_prefix in overlay_modes:
        out_fname = f"{file_prefix}_all_samples.tiff"
        render_sample_enrichment_overlay(
            df_sub=df_filt.copy(),
            out_path=OUT_DIR / out_fname,
            title_str=f"Distribution of Viruses – {mode_title} Overlay (All Samples)",
            n_total_cohort=N_TOTAL,
            strat_tag=SOURCE_FILE_FILTER,
            mode=mode,
        )
        out_fname_class = f"{file_prefix}_all_samples_by_class.tiff"
        render_specimen_class_overlay(
            df_cohort=df_filt.copy(),
            out_path=OUT_DIR / out_fname_class,
            title_str=f"Distribution of Viruses – {mode_title} by Specimen Class (All Samples)",
            n_total_cohort=N_TOTAL,
            strat_tag=SOURCE_FILE_FILTER,
            mode=mode,
        )


# ── Helper: demographics plot (vertical bars / histogram) ─────────────────────

def plot_demographics(data, col, level_label, filename, is_numeric=False):
    """
    Per-sample demographics: deduplicate to one row per sample_id first.
    Numeric columns (age, coverage) → histogram with subtle bar spacing.
    Categorical columns → vertical bar plot with square canvas (1800 x 1800 px) for collage matching.
    Top and right spines removed.
    """
    sample_df = data.drop_duplicates(subset="sample_id")

    # Exact square canvas (6.0 x 6.0 in, 1800 x 1800 px at 300 DPI) matching individual labeled plots
    fig, ax = plt.subplots(figsize=(6.0, 6.0), dpi=300)

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    if is_numeric:
        vals = pd.to_numeric(sample_df[col], errors="coerce").dropna()
        if vals.empty:
            plt.close(fig)
            print(f"  ⚠ no numeric data for {col}, skipping")
            return False
        if vals.nunique() <= 1:
            plt.close(fig)
            print(f"  ⚠ only one unique value for {col}, skipping")
            return False
        n_bins = min(20, max(8, int(np.sqrt(len(vals)))))
        ax.hist(vals, bins=n_bins, color=BLUE, edgecolor="none", rwidth=0.80)
        ax.set_xlabel(level_label, fontsize=10)
        ax.set_ylabel("Number of samples", fontsize=10)
        ax.tick_params(axis="x", labelsize=9, direction="out", length=3, width=0.6)
    else:
        cat_counts = sample_df[col].fillna("Unknown").value_counts()
        if len(cat_counts) <= 1:
            plt.close(fig)
            print(f"  ⚠ only one unique value for {col}, skipping")
            return False
        cat_counts = cat_counts.sort_values(ascending=False)
        x_labels = [str(x).replace("_", " ") for x in cat_counts.index]
        n_cats = len(cat_counts)

        # For 2-category plots (sex, phenotype, tissue): center bars with sleek width
        if n_cats == 2:
            x_pos = [-0.35, 0.35]
            bar_w = 0.45
            bars = ax.bar(x_pos, cat_counts.values,
                          color=BLUE, edgecolor="none", width=bar_w)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(x_labels, fontsize=9)
            ax.set_xlim(-1.0, 1.0)
        else:
            bar_w = 0.65
            bars = ax.bar(range(n_cats), cat_counts.values,
                          color=BLUE, edgecolor="none", width=bar_w)
            ax.set_xticks(range(n_cats))
            ax.set_xticklabels(x_labels, fontsize=9)

        max_val = cat_counts.max() if len(cat_counts) else 1
        for bar, cnt in zip(bars, cat_counts.values):
            h = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2, h + max_val * 0.015,
                f"{cnt}", ha="center", va="bottom", fontsize=8.5,
                color="black", fontweight="normal",
            )
        ax.set_ylabel("Number of samples", fontsize=10)
        ax.set_ylim(0, max_val * 1.15)
        ax.tick_params(axis="x", labelsize=9, direction="out", length=3, width=0.6)
        if any(len(lbl) > 8 for lbl in x_labels):
            plt.setp(ax.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
        else:
            plt.setp(ax.get_xticklabels(), rotation=0, ha="center")

    ax.tick_params(axis="y", labelsize=9, direction="out", length=3, width=0.6)
    ax.set_title(f"Demographics – {level_label}", fontsize=10.5)

    has_long_labels = not is_numeric and any(len(lbl) > 8 for lbl in x_labels)
    bot_margin = 0.22 if has_long_labels else 0.14
    fig.subplots_adjust(left=0.18, right=0.94, top=0.92, bottom=bot_margin)
    fig.text(0.5, 0.015, f"source: {SOURCE_FILE_FILTER}", ha="center",
             fontsize=8.0, fontstyle="normal", color="#333333")

    out = OUT_DIR / filename
    fig.savefig(out, dpi=300, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  ✓ saved {out}")
    return True


def plot_multi_virus_density(df_filt, out_path, title_str="Multi-Virus Copy Number Distribution", top_n=15, x_scale="log"):
    """
    Renders multi-virus smoothed copy-number density curve (resembling sketch):
    - Each virus is one smoothed line in a distinct color from PALETTE.
    - X-axis: copy number (log or linear).
    - Y-axis: counts / smoothed density.
    """
    if df_filt.empty:
        return

    # Filter to positive copy numbers
    df_pos = df_filt[df_filt["viral_copy_number"] > 0].copy()
    if df_pos.empty:
        return

    # Select top N viruses by positive sample count
    v_counts = df_pos.groupby("virus_accession")["sample_id"].nunique().sort_values(ascending=False)
    target_accessions = list(v_counts.index[:top_n]) if top_n > 0 else list(v_counts.index)

    fig, ax = plt.subplots(figsize=(9.0, 5.5), dpi=300)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.6)
    ax.spines["bottom"].set_linewidth(0.6)

    # Plot each virus
    for i, acc in enumerate(target_accessions):
        sub = df_pos[df_pos["virus_accession"] == acc]
        vals = sub["viral_copy_number"].values
        if len(vals) < 2:
            continue
        v_name = sub["virus_name_sanitized"].iloc[0] if "virus_name_sanitized" in sub.columns else acc
        color = PALETTE[i % len(PALETTE)]

        if x_scale == "log":
            vals_log = np.log10(vals[vals > 0])
            if len(vals_log) < 2:
                continue
            kde = gaussian_kde(vals_log, bw_method=0.4)
            x_grid = np.linspace(vals_log.min() - 0.2, vals_log.max() + 0.2, 200)
            y_dens = kde(x_grid) * len(vals)
            ax.plot(10**x_grid, y_dens, color=color, linewidth=2.0, label=f"{v_name} (n={len(vals)})")
        else:
            kde = gaussian_kde(vals, bw_method=0.4)
            x_grid = np.linspace(max(0, vals.min()), vals.max(), 200)
            y_dens = kde(x_grid) * len(vals)
            ax.plot(x_grid, y_dens, color=color, linewidth=2.0, label=f"{v_name} (n={len(vals)})")

    if x_scale == "log":
        ax.set_xscale("log")
        ax.set_xlabel("Viral Copy Number (Log Scale)", fontsize=10.5)
    else:
        ax.set_xlabel("Viral Copy Number", fontsize=10.5)

    ax.set_ylabel("Smoothed Positive Sample Count", fontsize=10.5)
    ax.set_title(title_str, fontsize=11.5, pad=10)
    ax.legend(bbox_to_anchor=(1.02, 0.5), loc="center left", frameon=False, fontsize=8.0)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  ✓ saved multi-virus density plot: {out_path}")


def plot_excluded_ubiquitous_viruses(df_filt, excluded_viruses, out_path, title_str="Excluded Ubiquitous Viruses Overview"):
    """
    Renders dedicated prevalence and copy number panel for excluded ubiquitous viruses.
    """
    if not excluded_viruses or df_filt.empty:
        return

    sub = df_filt[df_filt["virus_accession"].isin(excluded_viruses)].copy()
    if sub.empty:
        return

    n_v = len(excluded_viruses)
    fig, axes = plt.subplots(n_v, 2, figsize=(10.0, 3.0 * n_v), dpi=300, squeeze=False)

    for i, acc in enumerate(excluded_viruses):
        v_df = sub[sub["virus_accession"] == acc]
        v_name = v_df["virus_name_sanitized"].iloc[0] if not v_df.empty and "virus_name_sanitized" in v_df.columns else acc

        ax1 = axes[i, 0]
        ax2 = axes[i, 1]

        for a in [ax1, ax2]:
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
            a.spines["left"].set_linewidth(0.6)
            a.spines["bottom"].set_linewidth(0.6)

        # Panel 1: Reads distribution
        reads = v_df["virus_mapped_reads"].values if not v_df.empty else [0]
        ax1.hist(reads, bins=15, color=BLUE, rwidth=0.85)
        ax1.set_title(f"{v_name} ({acc}) – Read Counts", fontsize=9.5)
        ax1.set_xlabel("Mapped Reads", fontsize=8.5)
        ax1.set_ylabel("Number of Samples", fontsize=8.5)

        # Panel 2: Copy Number distribution
        cns = v_df[v_df["viral_copy_number"] > 0]["viral_copy_number"].values
        if len(cns) > 0:
            ax2.hist(cns, bins=15, color="#009E73", rwidth=0.85)
            ax2.set_xscale("log")
        ax2.set_title(f"{v_name} ({acc}) – Copy Number", fontsize=9.5)
        ax2.set_xlabel("Copy Number (Log Scale)", fontsize=8.5)
        ax2.set_ylabel("Number of Samples", fontsize=8.5)

    plt.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    print(f"  ✓ saved excluded ubiquitous viruses panel: {out_path}")

def parse_args():
    parser = argparse.ArgumentParser(description="Unified Viral Cohort Statistics & Visualizations (Module 5).")
    parser.add_argument("--input-report", required=True, help="Path to master viral report TSV file")
    parser.add_argument("--out-dir", default=".", help="Base output directory")
    parser.add_argument("--dataset", default="DATASET", help="Dataset identifier")
    parser.add_argument("--genome-build", default="hg38", help="Genome build")
    parser.add_argument("--target-phase", default="", help="Target phase filter (leave empty for auto)")
    parser.add_argument("--target-project", default="", help="Target project filter (leave empty for auto)")
    parser.add_argument("--strategies", nargs="*", default=["exogeneSR_viral_clean_filtered.sorted.flags.cram"], help="Target strategy CRAM names")
    parser.add_argument("--cohort-scope", default="both", choices=["both", "combined_all", "subsets_only", "target_only"], help="Cohort processing scope")
    parser.add_argument("--viral-short-names-file", default="HumanViral_Reference_02-07-2022_modified.renamed_map_short.tsv", help="Short name mapping TSV")
    parser.add_argument("--stats-filter-labels-file", default="", help="Optional labels file to filter clean cohort")
    parser.add_argument("--stats-generate-filtered-set", default="off", choices=["on", "off"], help="Generate clean vs all sample outputs")

    # Modular subsection switches
    parser.add_argument("--enable-summary-tables", default="on", choices=["on", "off"], help="5.1 Core Summary Tables")
    parser.add_argument("--enable-taxonomy-plots", default="on", choices=["on", "off"], help="5.2 Taxonomic Distribution Plots")
    parser.add_argument("--enable-broken-axis-taxonomy", default="on", choices=["on", "off"], help="Enable broken axes on taxonomy plots")
    parser.add_argument("--enable-family-background-overlay", default="on", choices=["on", "off"], help="Enable family overlay blocks")
    parser.add_argument("--sample-enrichment-plots", default="on", choices=["on", "off"], help="Enable sample enrichment overlay plots")
    parser.add_argument("--enable-demographics-plots", default="on", choices=["on", "off"], help="5.3 Demographics Plots")
    parser.add_argument("--enable-abundance-distributions", default="on", choices=["on", "off"], help="5.4 Abundance Distributions")
    parser.add_argument("--enable-multi-virus-density-plot", default="on", choices=["on", "off"], help="Enable multi-virus smoothed curve")
    parser.add_argument("--multi-density-top-n-viruses", type=int, default=15, help="Top N viruses for density plot")
    parser.add_argument("--multi-density-x-scale", default="log", choices=["log", "linear"], help="Density X scale")
    parser.add_argument("--reads-panel-a-loglog", default="on", choices=["on", "off"], help="Reads Panel A log-log")
    parser.add_argument("--reads-panel-b-log-y", default="on", choices=["on", "off"], help="Reads Panel B log-y")
    parser.add_argument("--copy-number-panel-a-loglog", default="on", choices=["on", "off"], help="Copy number Panel A log-log")
    parser.add_argument("--copy-number-panel-b-log-y", default="on", choices=["on", "off"], help="Copy number Panel B log-y")
    parser.add_argument("--log-scale-read-cutoff", type=int, default=30, help="Log scale read cutoff")
    parser.add_argument("--log-scale-copy-number-cutoff", type=float, default=0.1, help="Log scale copy number cutoff")
    parser.add_argument("--prelim-prevalence-cutoff-pct", type=float, default=5.0, help="Prevalence cutoff %")
    parser.add_argument("--prelim-mean-read-cutoff", type=float, default=6.0, help="Mean read cutoff")
    parser.add_argument("--enable-heatmaps", default="on", choices=["on", "off"], help="5.5 Cross-Sample Heatmaps")
    parser.add_argument("--heatmap-read-counts", default="off", choices=["on", "off"], help="Heatmap by reads")
    parser.add_argument("--heatmap-copy-number", default="off", choices=["on", "off"], help="Heatmap by copy number")
    parser.add_argument("--heatmap-min-reads-cutoff", type=int, default=3, help="Min reads for heatmap")
    parser.add_argument("--heatmap-stratify-by", default="", help="Heatmap stratification variable")
    parser.add_argument("--heatmap-exclude-viruses", nargs="*", default=[], help="Accessions to exclude from heatmaps")
    parser.add_argument("--enable-spatial-entropy-plots", default="on", choices=["on", "off"], help="5.6 Spatial Entropy Plots")
    parser.add_argument("--enable-coverage-qc-plots", default="on", choices=["on", "off"], help="5.7 Coverage Depth & Artifact Plots")
    parser.add_argument("--coverage-artifact-peak-fold-threshold", type=float, default=5.0, help="Peak threshold fold over median")
    parser.add_argument("--coverage-artifact-min-sample-recurrence", type=int, default=3, help="Min samples exhibiting peak to flag")
    parser.add_argument("--group-level", default="species", choices=["species", "genus", "family", "order", "realm", "accession", "none"], help="Taxonomic group level")
    parser.add_argument("--plots-dir", default="", help="Explicit plots dir")
    parser.add_argument("--stats-dir", default="", help="Explicit stats dir")
    return parser.parse_args()

def run_unified_stats(args):
    input_report_path = args.input_report
    base_out_dir = args.out_dir
    dataset = args.dataset or "DATASET"
    genome_build = args.genome_build or "hg38"
    stats_dir_name = args.stats_dir or f"{dataset}_{genome_build}_stats"
    plots_dir_name = args.plots_dir or f"{dataset}_{genome_build}_plots"

    if not os.path.exists(input_report_path):
        print(f"[ERROR] Input master report file not found: {input_report_path}")
        sys.exit(1)

    print(f"[START] Reading Master Viral Report: {input_report_path}")
    df_raw = pd.read_csv(input_report_path, sep="\t", low_memory=False, dtype=str)
    df_raw.columns = [c.lower().strip() for c in df_raw.columns]

    # Flexible column resolution
    def find_col(keywords, default_idx=0):
        for kw in keywords:
            for c in df_raw.columns:
                if kw in c:
                    return c
        return df_raw.columns[default_idx] if len(df_raw.columns) > default_idx else ""

    col_sample = find_col(["sample_id", "sample"], 0)
    col_accession = find_col(["virus_accession", "accession"], 1)
    col_reads = find_col(["virus_mapped_reads", "mapped_reads", "reads"], 3)
    col_cn = find_col(["viral_copy_number", "copy_number"], 8)
    col_name = find_col(["virus_name_sanitized", "virus_name", "name"], 9)
    col_source = find_col(["source_file", "strategy", "source"], 10)
    col_phase = find_col(["phase"], 11)
    col_project = find_col(["project"], 12)

    df = pd.DataFrame()
    df["sample_id"] = df_raw[col_sample].fillna("").astype(str).str.strip()
    df["virus_accession"] = df_raw[col_accession].fillna("").astype(str).str.strip()
    df["virus_mapped_reads"] = pd.to_numeric(df_raw[col_reads], errors="coerce").fillna(0).astype(int)
    df["viral_copy_number"] = pd.to_numeric(df_raw[col_cn], errors="coerce").fillna(0.0).astype(float)
    df["virus_name_sanitized"] = df_raw[col_name].fillna("").astype(str).str.strip()
    df["source_file"] = df_raw[col_source].fillna("unknown").astype(str).str.strip()
    df["phase"] = df_raw[col_phase].fillna("unknown").astype(str).str.strip()
    df["project"] = df_raw[col_project].fillna("base").astype(str).str.strip()
    df["project"] = df["project"].replace({"": "base", "none": "base", "0": "base"})

    for rank in ["species_taxid", "species_name", "genus_taxid", "genus_name", "family_taxid", "family_name", "order_taxid", "order_name", "realm_taxid", "realm_name"]:
        c = find_col([rank])
        df[rank] = df_raw[c].fillna("Unknown").astype(str).str.strip() if c else "Unknown"

    c_ent_h = find_col(["shannon_entropy", "entropy_h"])
    c_ent_j = find_col(["shannon_evenness", "evenness_j"])
    df["viral_shannon_entropy_h"] = pd.to_numeric(df_raw[c_ent_h], errors="coerce").fillna(0.0) if c_ent_h else 0.0
    df["viral_shannon_evenness_j"] = pd.to_numeric(df_raw[c_ent_j], errors="coerce").fillna(0.0) if c_ent_j else 0.0

    # Optional demographic columns
    for d_col in ["sex", "age", "tissue", "cell", "race", "clinic", "phenotype", "specimen"]:
        c = find_col([d_col])
        if c:
            df[d_col] = df_raw[c].fillna("Unknown").astype(str).str.strip()

    # Rescue taxonomy by keywords (papilloma, herpes, torque teno)
    mask_pap = df["virus_name_sanitized"].str.contains("papilloma", case=False, na=False)
    df.loc[mask_pap & (df["family_name"] == "Unknown"), "family_name"] = "Papillomaviridae"
    df.loc[mask_pap & (df["genus_name"] == "Unknown"), "genus_name"] = "Betapapillomavirus"
    df.loc[mask_pap & (df["order_name"] == "Unknown"), "order_name"] = "Zurhausenvirales"
    df.loc[mask_pap & (df["realm_name"] == "Unknown"), "realm_name"] = "Floreoviria"

    mask_herpes = df["virus_name_sanitized"].str.contains("herpes", case=False, na=False)
    df.loc[mask_herpes & (df["family_name"] == "Unknown"), "family_name"] = "Orthoherpesviridae"
    df.loc[mask_herpes & (df["order_name"] == "Unknown"), "order_name"] = "Herpesvirales"
    df.loc[mask_herpes & (df["realm_name"] == "Unknown"), "realm_name"] = "Duplodnaviria"

    mask_ttv = df["virus_name_sanitized"].str.contains("torque teno|anellovir", case=False, na=False)
    df.loc[mask_ttv & (df["family_name"] == "Unknown"), "family_name"] = "Anelloviridae"

    # Load short names map if provided or available
    short_map_dict = {}
    short_file_candidate = args.viral_short_names_file
    search_paths = [
        short_file_candidate,
        os.path.join(base_out_dir, "refs", short_file_candidate),
        os.path.join("/mnt/disks/staff/refs", short_file_candidate),
        os.path.join("refs", short_file_candidate),
        os.path.join("config/db_metadata", short_file_candidate),
        os.path.join("OLD", "HumanViral_Reference_02-07-2022_modified.renamed_map_short.tsv")
    ]
    resolved_short_map = None
    for p in search_paths:
        if p and os.path.exists(p):
            resolved_short_map = p
            break

    if resolved_short_map:
        try:
            sm_df = pd.read_csv(resolved_short_map, sep="\t")
            id_col = sm_df.columns[0]
            name_col = sm_df.columns[1] if len(sm_df.columns) > 1 else sm_df.columns[0]
            short_map_dict = dict(zip(sm_df[id_col].astype(str).str.strip(), sm_df[name_col].astype(str).str.strip()))
            print(f"[INFO] Loaded short names map from {resolved_short_map} ({len(short_map_dict):,} entries)")
        except Exception as e:
            print(f"[WARNING] Could not parse short names file {resolved_short_map}: {e}")

    if short_map_dict:
        df["virus_name_sanitized"] = df["virus_accession"].map(short_map_dict).fillna(df["virus_name_sanitized"])
    df["virus_name_sanitized"] = df["virus_name_sanitized"].apply(clean_virus_name)

    # Filter to requested strategy
    if args.strategies:
        df_strat = df[df["source_file"].isin(args.strategies)].copy()
        if df_strat.empty:
            print(f"[WARNING] No records match strategies {args.strategies}; using full report")
            df_strat = df.copy()
    else:
        df_strat = df.copy()

    # Determine execution subsets
    # If phase and/or project explicitly set, run ONLY that subset.
    # If empty, discover all distinct (phase, project) pairs + combined cohort.
    explicit_phase = args.target_phase.strip()
    explicit_project = args.target_project.strip()

    subsets = []
    if explicit_phase or explicit_project:
        subsets.append((explicit_phase, explicit_project, "target_only"))
    else:
        scope = args.cohort_scope.lower()
        if scope in ["combined_all", "both"]:
            subsets.append(("", "", "combined_all"))
        if scope in ["subsets_only", "both"]:
            pairs = df_strat.groupby(["phase", "project"]).size().index.tolist()
            for ph, pr in pairs:
                if ph and ph.lower() != "unknown":
                    subsets.append((ph, pr, "subset"))

    print(f"[INFO] Configured execution subsets to process: {len(subsets)}")

    # Process each subset
    for cur_phase, cur_proj, subset_type in subsets:
        if subset_type == "combined_all":
            df_cur = df_strat.copy()
            cur_phase_tag = "all"
            cur_proj_tag = "all"
            sub_prefix = f"{dataset}_{genome_build}_all"
            cur_stats_dir = os.path.join(base_out_dir, stats_dir_name)
            cur_plots_dir = os.path.join(base_out_dir, plots_dir_name)
            cov_plots_dir = os.path.join(cur_plots_dir, "coverage_plots")
        else:
            # Normalize phase naming (phase1 vs 1)
            clean_p = cur_phase.replace("phase", "").lower()
            df_cur = df_strat[
                (df_strat["phase"].str.lower().str.replace("phase", "") == clean_p) &
                (df_strat["project"].str.lower() == cur_proj.lower())
            ].copy()
            if clean_p.isdigit():
                cur_phase_tag = f"phase{clean_p}"
            else:
                cur_phase_tag = cur_phase

            cur_proj_tag = cur_proj if cur_proj else "base"
            sub_prefix = f"{dataset}_{genome_build}_{cur_phase_tag}_{cur_proj_tag}"
            cur_stats_dir = os.path.join(base_out_dir, stats_dir_name, cur_phase_tag, cur_proj_tag)
            cur_plots_dir = os.path.join(base_out_dir, plots_dir_name, cur_phase_tag, cur_proj_tag)
            cov_plots_dir = os.path.join(cur_plots_dir, "coverage_plots")

        os.makedirs(cur_stats_dir, exist_ok=True)
        os.makedirs(cur_plots_dir, exist_ok=True)
        os.makedirs(cov_plots_dir, exist_ok=True)

        print(f"\n==================================================================")
        print(f"Processing Scope: Phase='{cur_phase_tag}', Project='{cur_proj_tag}'")
        print(f"Stats Directory:  {cur_stats_dir}")
        print(f"Plots Directory:  {cur_plots_dir}")
        print(f"Rows in subset:   {len(df_cur):,}")
        print(f"==================================================================")

        if df_cur.empty:
            print("[INFO] Subset is empty; continuing.")
            continue

        n_total_cohort = df_cur["sample_id"].nunique()

        # Filter out rows with no virus detected for viral-level statistics
        df_viral = df_cur[df_cur["virus_accession"].notna() & (df_cur["virus_accession"] != "") & (df_cur["virus_accession"].str.lower() != "none") & (df_cur["virus_mapped_reads"] > 0)].copy()

        # ----------------------------------------------------------------------
        # 5.1: Core Summary Tables
        # ----------------------------------------------------------------------
        if args.enable_summary_tables == "on":
            print("── [5.1] Generating Core Summary Tables ──")
            tsv_summary = os.path.join(cur_stats_dir, f"{sub_prefix}_stats_summary.tsv")
            md_summary = os.path.join(cur_stats_dir, f"{sub_prefix}_stats_summary.md")
            
            # Write master metrics
            with open(tsv_summary, "w") as f_tsv, open(md_summary, "w") as f_md:
                header = ["Metric", "Value"]
                f_tsv.write("\t".join(header) + "\n")
                f_md.write(f"# Viral QC & Cohort Summary: {sub_prefix}\n\n")
                f_md.write("| Metric | Value |\n| :--- | :--- |\n")

                metrics = [
                    ("Dataset", dataset),
                    ("Genome Build", genome_build),
                    ("Phase", cur_phase_tag),
                    ("Project", cur_proj_tag),
                    ("Total Cohort Samples", str(n_total_cohort)),
                    ("Viral Positive Samples", str(df_viral["sample_id"].nunique())),
                    ("Cohort Viral Prevalence (%)", f"{df_viral['sample_id'].nunique() / float(n_total_cohort) * 100.0:.2f}%" if n_total_cohort > 0 else "0%"),
                    ("Total Filtered Viral Reads", f"{df_viral['virus_mapped_reads'].sum():,}"),
                    ("Unique Viral Taxa Detected", str(df_viral["virus_accession"].nunique())),
                ]
                for m, val in metrics:
                    f_tsv.write(f"{m}\t{val}\n")
                    f_md.write(f"| {m} | {val} |\n")
            print(f"  ✓ Written summary TSV & MD: {tsv_summary}")

        # ----------------------------------------------------------------------
        # 5.2: Taxonomic Distribution & Landscape Plots
        # ----------------------------------------------------------------------
        if args.enable_taxonomy_plots == "on" and not df_viral.empty:
            print("── [5.2] Generating Taxonomic Distribution Plots ──")
            enable_broken = (args.enable_broken_axis_taxonomy == "on")
            tax_levels = [
                ("virus_accession",  "virus_name_sanitized", "Accession",  f"{sub_prefix}_taxonomy_accession"),
                ("species_taxid",    "species_name",         "Species",    f"{sub_prefix}_taxonomy_species"),
                ("genus_taxid",      "genus_name",           "Genus",      f"{sub_prefix}_taxonomy_genus"),
                ("family_taxid",     "family_name",          "Family",     f"{sub_prefix}_taxonomy_family"),
                ("order_taxid",      "order_name",           "Order",      f"{sub_prefix}_taxonomy_order"),
                ("realm_taxid",      "realm_name",           "Realm",      f"{sub_prefix}_taxonomy_realm"),
            ]
            for taxid_col, name_col, lvl_label, fname_base in tax_levels:
                plot_taxonomy_level(
                    data=df_viral,
                    taxid_col=taxid_col,
                    name_col=name_col,
                    level_label=lvl_label,
                    save_path_prefix=os.path.join(cur_plots_dir, fname_base),
                    n_total_cohort=n_total_cohort,
                    strat_tag=args.strategies[0] if args.strategies else "clean_flags",
                    enable_broken=enable_broken
                )

            if args.enable_family_background_overlay == "on":
                plot_accession_family_overlay(
                    data=df_viral,
                    save_path_prefix=os.path.join(cur_plots_dir, f"{sub_prefix}_taxonomy_accession_family_overlay"),
                    n_total_cohort=n_total_cohort,
                    strat_tag=args.strategies[0] if args.strategies else "clean_flags",
                    enable_broken=enable_broken
                )

            if args.sample_enrichment_plots == "on":
                overlay_modes = [
                    ("samples", "Sample Enrichment", f"{sub_prefix}_sample_enrichment_overlay_samples"),
                    ("reads",   "Read Abundance",    f"{sub_prefix}_sample_enrichment_overlay_reads"),
                    ("copynumber", "Viral Copy Number", f"{sub_prefix}_sample_enrichment_overlay_copynumber"),
                ]
                for mode, mode_title, pfx in overlay_modes:
                    render_sample_enrichment_overlay(
                        df_sub=df_viral.copy(),
                        out_path=pathlib.Path(cur_plots_dir) / f"{pfx}_all_samples.tiff",
                        title_str=f"Distribution of Viruses – {mode_title} Overlay (All Samples)",
                        n_total_cohort=n_total_cohort,
                        strat_tag=args.strategies[0] if args.strategies else "clean_flags",
                        mode=mode
                    )

        # ----------------------------------------------------------------------
        # 5.3: Cohort Demographics Plots (Strictly Non-Broken Axes)
        # ----------------------------------------------------------------------
        if args.enable_demographics_plots == "on":
            print("── [5.3] Generating Cohort Demographics Plots ──")
            demo_cols = [
                ("sex",              "Sex",       f"{sub_prefix}_demo_sex.tiff",       False),
                ("age",              "Age",       f"{sub_prefix}_demo_age.tiff",       True),
                ("tissue",           "Tissue",    f"{sub_prefix}_demo_tissue.tiff",    False),
                ("cell",             "Cell Type", f"{sub_prefix}_demo_cell.tiff",      False),
                ("race",             "Race",      f"{sub_prefix}_demo_race.tiff",      False),
                ("phenotype",        "Phenotype", f"{sub_prefix}_demo_phenotype.tiff", False),
                ("specimen",         "Specimen",  f"{sub_prefix}_demo_specimen.tiff",  False),
            ]
            for col, lbl, fn, is_num in demo_cols:
                if col in df_cur.columns:
                    plot_demographics(df_cur, col, lbl, pathlib.Path(cur_plots_dir) / fn, is_numeric=is_num)

        # ----------------------------------------------------------------------
        # 5.4: Abundance & Multi-Virus Density Distributions
        # ----------------------------------------------------------------------
        if args.enable_abundance_distributions == "on" and not df_viral.empty:
            print("── [5.4] Generating Abundance & Multi-Virus Density Distributions ──")
            if args.enable_multi_virus_density_plot == "on":
                dens_path = os.path.join(cur_plots_dir, f"{sub_prefix}_copy_number_multi_density.tiff")
                plot_multi_virus_density(
                    df_filt=df_viral,
                    out_path=dens_path,
                    title_str=f"Multi-Virus Copy Number Distribution ({cur_phase_tag} {cur_proj_tag})",
                    top_n=args.multi_density_top_n_viruses,
                    x_scale=args.multi_density_x_scale
                )

        # ----------------------------------------------------------------------
        # 5.5: Cross-Sample Heatmaps & Ubiquitous Virus Reporting
        # ----------------------------------------------------------------------
        if args.enable_heatmaps == "on":
            print("── [5.5] Generating Heatmaps & Ubiquitous Virus Distribution ──")
            # If ubiquitous viruses excluded, render dedicated panel
            if args.heatmap_exclude_viruses:
                ex_path = os.path.join(cur_plots_dir, f"{sub_prefix}_excluded_ubiquitous_viruses_distribution.tiff")
                plot_excluded_ubiquitous_viruses(
                    df_filt=df_viral,
                    excluded_viruses=args.heatmap_exclude_viruses,
                    out_path=ex_path,
                    title_str=f"Excluded Ubiquitous Viruses Distribution ({cur_phase_tag} {cur_proj_tag})"
                )

        # ----------------------------------------------------------------------
        # 5.6: Spatial Uniformity & Entropy Plots
        # ----------------------------------------------------------------------
        if args.enable_spatial_entropy_plots == "on" and not df_viral.empty:
            print("── [5.6] Generating Spatial Entropy & Evenness Boxplots ──")
            df_ent = df_viral[(df_viral["viral_shannon_entropy_h"] > 0) | (df_viral["viral_shannon_evenness_j"] > 0)]
            if not df_ent.empty:
                # 1. Multi-virus comparative boxplot: Entropy H
                fig, ax = plt.subplots(figsize=(10, max(5.0, df_ent["virus_name_sanitized"].nunique() * 0.35)), dpi=300)
                v_groups_h = df_ent.groupby("virus_name_sanitized")["viral_shannon_entropy_h"].apply(list)
                names_h = list(v_groups_h.index)
                data_h = [v_groups_h[n] for n in names_h]
                ax.boxplot(data_h, vert=False, labels=names_h, patch_artist=True, boxprops=dict(facecolor="#4472C4", alpha=0.7))
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.set_xlabel("Shannon Entropy (H)", fontsize=10.5)
                ax.set_title(f"Viral Spatial Shannon Entropy (H) – {cur_phase_tag} {cur_proj_tag}", fontsize=11.5)
                plt.tight_layout()
                h_path = os.path.join(cur_plots_dir, f"{sub_prefix}_shannon_entropy_H.tiff")
                fig.savefig(h_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
                plt.close(fig)
                print(f"  ✓ saved {h_path}")

                # 2. Multi-virus comparative boxplot: Evenness J
                fig, ax = plt.subplots(figsize=(10, max(5.0, df_ent["virus_name_sanitized"].nunique() * 0.35)), dpi=300)
                v_groups_j = df_ent.groupby("virus_name_sanitized")["viral_shannon_evenness_j"].apply(list)
                names_j = list(v_groups_j.index)
                data_j = [v_groups_j[n] for n in names_j]
                ax.boxplot(data_j, vert=False, labels=names_j, patch_artist=True, boxprops=dict(facecolor="#ED7D31", alpha=0.7))
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.set_xlim(0.0, 1.05)
                ax.set_xlabel("Pielou's Evenness (J)", fontsize=10.5)
                ax.set_title(f"Viral Spatial Pielou's Evenness (J) – {cur_phase_tag} {cur_proj_tag}", fontsize=11.5)
                plt.tight_layout()
                j_path = os.path.join(cur_plots_dir, f"{sub_prefix}_pielou_evenness_J.tiff")
                fig.savefig(j_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
                plt.close(fig)
                print(f"  ✓ saved {j_path}")

        # ----------------------------------------------------------------------
        # 5.7: Coverage Depth, Hotspots & Artifact Exclusion
        # ----------------------------------------------------------------------
        if args.enable_coverage_qc_plots == "on" and not df_viral.empty:
            print("── [5.7] Generating Pooled Cohort Coverage & Artifact Tables ──")
            # Discover and aggregate per-hit binned coverage files
            # Check candidate report directories:
            report_glob_base = os.path.join(base_out_dir, f"{dataset}_{genome_build}_reports")
            artifact_rows = []
            
            # Analyze each detected virus accession in this subset
            unique_accs = df_viral["virus_accession"].unique()
            for acc in unique_accs:
                # Look for {SAMPLE}_{acc}_genome_bins_coverage.tsv
                pattern = os.path.join(report_glob_base, "**", f"*_{acc}_genome_bins_coverage.tsv")
                bin_files = glob.glob(pattern, recursive=True)

                if not bin_files:
                    continue

                bin_data_list = []
                for bf in bin_files:
                    try:
                        bdf = pd.read_csv(bf, sep="\t")
                        if "bin_index" in bdf.columns and "mean_depth" in bdf.columns:
                            bin_data_list.append(bdf)
                    except Exception:
                        continue

                if not bin_data_list:
                    continue

                all_bins_df = pd.concat(bin_data_list, ignore_index=True)
                pooled_summary = all_bins_df.groupby("bin_index").agg(
                    bin_start=("bin_start", "min"),
                    bin_end=("bin_end", "max"),
                    mean_depth=("mean_depth", "mean"),
                    median_depth=("mean_depth", "median"),
                    sample_count=("mean_depth", "count")
                ).reset_index()

                # Write pooled coverage TSV
                pooled_tsv_path = os.path.join(cur_stats_dir, f"{sub_prefix}_{acc}_pooled_genome_coverage.tsv")
                pooled_summary.to_csv(pooled_tsv_path, sep="\t", index=False)

                # Render pooled coverage plot into coverage_plots/
                pooled_plot_path = os.path.join(cov_plots_dir, f"{sub_prefix}_{acc}_pooled_genome_coverage.tiff")
                fig, ax = plt.subplots(figsize=(8.5, 4.5), dpi=300)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                bars = ax.bar(pooled_summary["bin_index"], pooled_summary["mean_depth"], color="#0072B2", edgecolor="none", width=0.85)
                ax.set_xlabel("Genomic Percentile Bin (1 - 10)", fontsize=10)
                ax.set_ylabel("Mean Pooled Read Depth", fontsize=10)
                v_label = df_viral[df_viral["virus_accession"] == acc]["virus_name_sanitized"].iloc[0] if not df_viral[df_viral["virus_accession"] == acc].empty else acc
                ax.set_title(f"Pooled Viral Genome Coverage: {v_label} ({acc})\nPhase: {cur_phase_tag} | Project: {cur_proj_tag}", fontsize=11)
                plt.tight_layout()
                fig.savefig(pooled_plot_path, dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
                plt.close(fig)

                # Check for recurrent artifact peaks
                med_overall = pooled_summary["mean_depth"].median()
                threshold = med_overall * args.coverage_artifact_peak_fold_threshold
                for row_b in pooled_summary.itertuples():
                    if row_b.mean_depth >= threshold and row_b.sample_count >= args.coverage_artifact_min_sample_recurrence:
                        artifact_rows.append([
                            dataset, genome_build, cur_phase_tag, cur_proj_tag, acc, v_label,
                            str(row_b.bin_index), str(row_b.bin_start), str(row_b.bin_end),
                            f"{row_b.mean_depth:.2f}", f"{med_overall:.2f}", str(row_b.sample_count),
                            "FLAGGED_HOTSPOT_ARTIFACT"
                        ])

            # Output cohort coverage artifact summary table
            artifact_tsv_path = os.path.join(cur_stats_dir, f"{sub_prefix}_cohort_coverage_artifact_regions.tsv")
            with open(artifact_tsv_path, "w") as f_art:
                f_art.write("\t".join([
                    "dataset", "genome_build", "phase", "project", "virus_accession", "virus_name",
                    "bin_index", "bin_start", "bin_end", "bin_mean_depth", "genome_median_depth",
                    "recurrence_sample_count", "artifact_status"
                ]) + "\n")
                for r in artifact_rows:
                    f_art.write("\t".join(r) + "\n")
            print(f"  ✓ Written artifact summary table: {artifact_tsv_path} ({len(artifact_rows)} regions flagged)")


def main():
    args = parse_args()
    run_unified_stats(args)
    from datetime import datetime
    print(f"\n[COMPLETED] Execution finished at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
