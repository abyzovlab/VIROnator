#!/usr/bin/env python3
"""
Flag Comparison Script for VIROnator
Compares read pair counts across 3 SAM flag filtering methods for the clean strategy:
1. original_command_flags_reads          : Original command WITH -f 2 flag
2. original_command_flags_modified_reads : Original command WITHOUT -f 2 flag
3. manual_flags_reads                    : Manual gawk bitwise flag check

Handles multi-phase and multi-project sample occurrences automatically.
"""

import os
import sys
import glob
import argparse
import subprocess


def parse_args():
    parser = argparse.ArgumentParser(description="Systematic SAM Flag Comparison for VIROnator Clean Strategy")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--phase", required=True, help="Phase number or tag")
    parser.add_argument("--project", default="base", help="Project name (default: base)")
    parser.add_argument("--vironator-dir", required=True, help="Path to sample VIROnator output directory")
    parser.add_argument("--viral-bed", required=True, help="Path to viral target regions BED file")
    parser.add_argument("--out-file", required=True, help="Output TSV report file path")
    return parser.parse_args()


def run_samtools_count(cmd):
    """Executes a samtools pipeline command and returns integer read count."""
    try:
        res = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        val = res.decode("utf-8").strip()
        return int(val) if val else 0
    except Exception:
        return 0


def extract_phase_project_from_path(cram_path):
    """
    Parses phase and project from path structure:
    .../SSC_hg38_vironator/phase1/Project_A/SS0012978/exogeneSR_viral_clean.cram
    """
    norm_path = os.path.normpath(cram_path)
    parts = norm_path.split(os.sep)
    
    phase = "unknown"
    project = "base"
    
    for i, p in enumerate(parts):
        if p.startswith("phase"):
            phase = p
            if i + 2 < len(parts) - 1: # Has project folder
                project = parts[i + 1]
            break
            
    return phase, project


def find_all_clean_crams(sample_id, vironator_dir):
    """
    Locates all clean CRAM files for sample_id.
    If vironator_dir contains the file directly, returns [ (phase, project, path) ].
    If not, searches up parent directories to find all phase/project instances of sample_id.
    """
    results = []
    candidates = [
        "exogeneSR_viral_clean.cram",
        "exogeneSR_viral_clean_filtered.sorted.cram",
        "exogeneSR_viral_clean_filtered.cram"
    ]
    
    # Check direct directory first
    if os.path.exists(vironator_dir):
        for c in candidates:
            fp = os.path.join(vironator_dir, c)
            if os.path.exists(fp):
                results.append(fp)
                break
                
    if not results and os.path.exists(vironator_dir):
        for f in os.listdir(vironator_dir):
            if f.endswith(".cram") and "clean" in f and not "flags" in f:
                results.append(os.path.join(vironator_dir, f))
                break

    # If nothing found directly, search parent vironator directory tree for sample_id
    if not results:
        curr = os.path.abspath(vironator_dir)
        # Walk up to find top-level vironator directory (containing 'phase' folders)
        top_dir = None
        while curr and curr != os.path.dirname(curr):
            if "vironator" in os.path.basename(curr).lower():
                top_dir = curr
                break
            curr = os.path.dirname(curr)
            
        if top_dir and os.path.exists(top_dir):
            pattern = os.path.join(top_dir, "**", sample_id, "*.cram")
            for fp in glob.glob(pattern, recursive=True):
                fname = os.path.basename(fp)
                if "clean" in fname and not "flags" in fname:
                    results.append(fp)

    return results


def process_cram(in_cram, viral_bed):
    """Runs 3 samtools commands on input clean CRAM file."""
    # Command 1: Original WITH -f 2
    cmd_orig = (
        f"samtools view -F 2048 -F 12 -f 2 -q 20 \"{in_cram}\" --region-file \"{viral_bed}\" "
        f"| awk '($9 >= 150 && $9 <= 1500) || ($9 <= -150 && $9 >= -1500)' "
        f"| awk '$6==\"151M\" || $6==\"150M\"' "
        f"| awk '$7==\"=\"' "
        f"| cut -f1 | sort | uniq -c | awk '$1==\"2\"' | wc -l"
    )
    
    # Command 2: Original WITHOUT -f 2
    cmd_modif = (
        f"samtools view -F 2048 -F 12 -q 20 \"{in_cram}\" --region-file \"{viral_bed}\" "
        f"| awk '($9 >= 150 && $9 <= 1500) || ($9 <= -150 && $9 >= -1500)' "
        f"| awk '$6==\"151M\" || $6==\"150M\"' "
        f"| awk '$7==\"=\"' "
        f"| cut -f1 | sort | uniq -c | awk '$1==\"2\"' | wc -l"
    )
    
    # Command 3: Manual gawk bitwise flag check
    cmd_manual = (
        f"samtools view -F 2048 -F 12 -q 20 \"{in_cram}\" --region-file \"{viral_bed}\" "
        f"| gawk '((and($2,97) && !and($2,16) && ($9 >= 150 && $9 <= 1500)) || (and($2,81) && !and($2,32) && ($9 >= -1500 && $9 <= -150)) || (and($2,145) && !and($2,32) && ($9 >= 150 && $9 <= 1500)) || (and($2,161) && !and($2,16) && ($9 >= -1500 && $9 <= -150)))' "
        f"| awk '$6==\"151M\" || $6==\"150M\"' "
        f"| awk '$7==\"=\"' "
        f"| cut -f1 | sort | uniq -c | awk '$1==\"2\"' | wc -l"
    )
    
    n_orig = run_samtools_count(cmd_orig)
    n_modif = run_samtools_count(cmd_modif)
    n_manual = run_samtools_count(cmd_manual)
    
    return n_orig, n_modif, n_manual


def main():
    args = parse_args()
    
    fallback_phase = f"phase{args.phase}" if not str(args.phase).startswith("phase") else str(args.phase)
    fallback_project = str(args.project) if args.project else "base"
    sample_id = str(args.sample_id)
    strategy_str = "clean"
    
    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    header = [
        "phase",
        "project",
        "sample",
        "strategy",
        "original_command_flags_reads",
        "original_command_flags_modified_reads",
        "manual_flags_reads"
    ]
    
    clean_crams = find_all_clean_crams(sample_id, args.vironator_dir)
    rows = []

    if not clean_crams or not os.path.exists(args.viral_bed):
        rows.append([fallback_phase, fallback_project, sample_id, strategy_str, "0", "0", "0"])
    else:
        for in_cram in clean_crams:
            parsed_phase, parsed_proj = extract_phase_project_from_path(in_cram)
            phase_val = parsed_phase if parsed_phase != "unknown" else fallback_phase
            proj_val = parsed_proj if parsed_proj else fallback_project
            
            n_orig, n_modif, n_manual = process_cram(in_cram, args.viral_bed)
            rows.append([phase_val, proj_val, sample_id, strategy_str, str(n_orig), str(n_modif), str(n_manual)])

    with open(args.out_file, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Flag comparison completed for sample {sample_id} ({len(rows)} instances processed).")


if __name__ == "__main__":
    main()
