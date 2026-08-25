#!/usr/bin/env python3
"""
Flag Comparison Script for VIROnator
Compares read pair counts across 3 SAM flag filtering methods for the clean strategy:
1. original_command_flags_reads          : Original command WITH -f 2 flag
2. original_command_flags_modified_reads : Original command WITHOUT -f 2 flag
3. manual_flags_reads                    : Manual gawk bitwise flag check
"""

import os
import sys
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


def find_clean_cram(vironator_dir):
    """Locates the input CRAM file for the clean strategy in vironator_dir."""
    candidates = [
        "exogeneSR_viral_clean.cram",
        "exogeneSR_viral_clean_filtered.sorted.cram",
        "exogeneSR_viral_clean_filtered.cram"
    ]
    for c in candidates:
        full_path = os.path.join(vironator_dir, c)
        if os.path.exists(full_path):
            return full_path
    
    # Fallback search for any CRAM with 'clean' in its filename
    if os.path.exists(vironator_dir):
        for f in os.listdir(vironator_dir):
            if f.endswith(".cram") and "clean" in f and not "flags" in f:
                return os.path.join(vironator_dir, f)
    return None


def main():
    args = parse_args()
    
    phase_str = str(args.phase)
    proj_str = str(args.project) if args.project else "base"
    sample_id = str(args.sample_id)
    strategy_str = "clean"
    
    out_dir = os.path.dirname(args.out_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    in_cram = find_clean_cram(args.vironator_dir)
    
    if not in_cram or not os.path.exists(in_cram):
        print(f"[WARNING] Input clean CRAM not found in {args.vironator_dir}. Writing zeros.")
        n_orig, n_modif, n_manual = 0, 0, 0
    elif not os.path.exists(args.viral_bed):
        print(f"[WARNING] Viral BED file not found: {args.viral_bed}. Writing zeros.")
        n_orig, n_modif, n_manual = 0, 0, 0
    else:
        # Command 1: Original WITH -f 2
        cmd_orig = (
            f"samtools view -F 2048 -F 12 -f 2 -q 20 \"{in_cram}\" --region-file \"{args.viral_bed}\" "
            f"| awk '($9 >= 150 && $9 <= 1500) || ($9 <= -150 && $9 >= -1500)' "
            f"| awk '$6==\"151M\" || $6==\"150M\"' "
            f"| awk '$7==\"=\"' "
            f"| cut -f1 | sort | uniq -c | awk '$1==\"2\"' | wc -l"
        )
        
        # Command 2: Original WITHOUT -f 2
        cmd_modif = (
            f"samtools view -F 2048 -F 12 -q 20 \"{in_cram}\" --region-file \"{args.viral_bed}\" "
            f"| awk '($9 >= 150 && $9 <= 1500) || ($9 <= -150 && $9 >= -1500)' "
            f"| awk '$6==\"151M\" || $6==\"150M\"' "
            f"| awk '$7==\"=\"' "
            f"| cut -f1 | sort | uniq -c | awk '$1==\"2\"' | wc -l"
        )
        
        # Command 3: Manual gawk bitwise flag check
        cmd_manual = (
            f"samtools view -F 2048 -F 12 -q 20 \"{in_cram}\" --region-file \"{args.viral_bed}\" "
            f"| gawk '((and($2,97) && !and($2,16) && ($9 >= 150 && $9 <= 1500)) || (and($2,81) && !and($2,32) && ($9 >= -1500 && $9 <= -150)) || (and($2,145) && !and($2,32) && ($9 >= 150 && $9 <= 1500)) || (and($2,161) && !and($2,16) && ($9 >= -1500 && $9 <= -150)))' "
            f"| awk '$6==\"151M\" || $6==\"150M\"' "
            f"| awk '$7==\"=\"' "
            f"| cut -f1 | sort | uniq -c | awk '$1==\"2\"' | wc -l"
        )
        
        n_orig = run_samtools_count(cmd_orig)
        n_modif = run_samtools_count(cmd_modif)
        n_manual = run_samtools_count(cmd_manual)

    header = [
        "phase",
        "project",
        "sample",
        "strategy",
        "original_command_flags_reads",
        "original_command_flags_modified_reads",
        "manual_flags_reads"
    ]
    
    row = [
        phase_str,
        proj_str,
        sample_id,
        strategy_str,
        str(n_orig),
        str(n_modif),
        str(n_manual)
    ]

    with open(args.out_file, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(row) + "\n")

    print(f"[SUCCESS] Flag comparison completed for {sample_id}: WITH_F2={n_orig}, WITHOUT_F2={n_modif}, MANUAL_FLAGS={n_manual}")


if __name__ == "__main__":
    main()
