#!/usr/bin/env python3
"""
build_refinement_ref.py
Combined Reference Builder & Bowtie2 Indexer for VIROnator Refinement Module

1. Concatenates Human (GRCh38), Mouse, Plasmids, and downloaded RefSeq viral FASTAs.
2. Writes date-stamped combined reference to /mnt/disks/staff/refs/combined_human_mouse_plasmids_refinement_viral_<DATE>.fa.
3. Checks if index already exists (skips re-indexing if present).
4. Generates samtools faidx (.fai) and bowtie2-build (.bt2 / .bt2l) indices.
"""

import os
import sys
import glob
import argparse
import subprocess
import time


def parse_args():
    parser = argparse.ArgumentParser(description="Build date-stamped combined reference and Bowtie2 index for refinement module.")
    parser.add_argument("--ref-human", required=True, help="Path to Human reference FASTA (GRCh38)")
    parser.add_argument("--ref-mouse", default="", help="Path to Mouse reference FASTA")
    parser.add_argument("--ref-plasmids", default="", help="Path to Plasmids reference FASTA")
    parser.add_argument("--download-dir", default="/mnt/disks/staff/refs/ncbi_download", help="Directory containing downloaded RefSeq viral FASTAs")
    parser.add_argument("--output-ref", default="/mnt/disks/staff/refs/combined_human_mouse_plasmids_refinement_viral_2026-08-29.fa", help="Path to date-stamped output combined reference FASTA")
    parser.add_argument("--bowtie2-build-bin", default="bowtie2-build", help="Path to bowtie2-build binary")
    return parser.parse_args()


def main():
    args = parse_args()

    out_fasta = args.output_ref
    out_dir = os.path.dirname(out_fasta)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Check if combined FASTA and Bowtie2 index already exist
    bt2_idx1 = out_fasta + ".1.bt2"
    bt2l_idx1 = out_fasta + ".1.bt2l"
    fai_idx = out_fasta + ".fai"

    if os.path.exists(out_fasta) and os.path.getsize(out_fasta) > 0 and (os.path.exists(bt2_idx1) or os.path.exists(bt2l_idx1)):
        print(f"[CACHE HIT] Date-stamped combined reference and Bowtie2 index already exist: {out_fasta}")
        print("[CACHE HIT] Skipping reference concatenation and indexing.")
        return

    print(f"[BUILD] Assembling combined reference FASTA: {out_fasta} ...", flush=True)

    # Gather viral fastas from download_dir (prefer .renamed.fa sanitized files)
    viral_fastas = []
    if os.path.exists(args.download_dir):
        renamed_files = glob.glob(os.path.join(args.download_dir, "*.renamed.fa"))
        if renamed_files:
            viral_fastas = renamed_files
        else:
            viral_fastas = [f for f in glob.glob(os.path.join(args.download_dir, "*.fa")) if not f.endswith(".renamed.fa")]

    print(f"[INFO] Including {len(viral_fastas)} downloaded RefSeq viral genome FASTAs.")

    # Concatenate Human + Mouse + Plasmids + Downloaded Viral FASTAs
    source_files = []
    if os.path.exists(args.ref_human):
        source_files.append(args.ref_human)
    if args.ref_mouse and os.path.exists(args.ref_mouse):
        source_files.append(args.ref_mouse)
    if args.ref_plasmids and os.path.exists(args.ref_plasmids):
        source_files.append(args.ref_plasmids)
    source_files.extend(viral_fastas)

    with open(out_fasta, "w", encoding="utf-8") as out_f:
        for sf in source_files:
            print(f"  Appending {os.path.basename(sf)} ...", flush=True)
            with open(sf, "r", encoding="utf-8", errors="replace") as in_f:
                for line in in_f:
                    out_f.write(line)

    print(f"[SUCCESS] Combined reference FASTA written: {out_fasta} ({os.path.getsize(out_fasta)} bytes)")

    # samtools faidx
    if not os.path.exists(fai_idx):
        print(f"[INDEX] Running samtools faidx {out_fasta} ...", flush=True)
        subprocess.run(f"samtools faidx \"{out_fasta}\"", shell=True, check=True)

    # bowtie2-build
    if not (os.path.exists(bt2_idx1) or os.path.exists(bt2l_idx1)):
        print(f"[INDEX] Running bowtie2-build on {out_fasta} (this may take a few minutes) ...", flush=True)
        bt2_cmd = f"\"{args.bowtie2-build-bin}\" --threads 16 \"{out_fasta}\" \"{out_fasta}\""
        subprocess.run(bt2_cmd, shell=True, check=True)

    print(f"[SUCCESS] Combined reference FASTA and Bowtie2 index construction complete: {out_fasta}")


if __name__ == "__main__":
    main()
