#!/usr/bin/env python3
"""
standalone_igv_snapshot_generator.py

================================================================================
INSTRUCTIONS & HOW TO RUN (STANDALONE IGV SNAPSHOT GENERATOR)
================================================================================

This script generates an IGV batch script (.igv) and automatically renders PNG 
snapshots for all sample-virus pairs with non-zero read differences (difference != 0).

3 CRAM tracks are loaded simultaneously per sample-virus hit:
  1. flags      : exogeneSR_viral_clean_filtered.sorted.flags.cram
  2. noflags    : exogeneSR_viral_clean_filtered.sorted.cram
  3. additional : exogeneSR_viral_clean_filtered.sorted.flags.additional.cram

Sets locus to full virus genome length (goto <virus_accession>) without zooming in.

--------------------------------------------------------------------------------
ONE-LINER COMMAND:
--------------------------------------------------------------------------------

python3 scripts/standalone_igv_snapshot_generator.py --diff-tsv MCBiobank_hg38_comparisons/MCBiobank_hg38_master_report_cleans_common_flags_vs_noflags_diff.tsv --vironator-dir ./CRAMs --ref-genome ./REFS/HumanViral_Reference_02-07-2022_modified.renamed.fa --igv-binary igv.sh --run-igv

================================================================================
"""

import argparse
import os
import sys
import subprocess
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Standalone IGV batch script generator and PNG snapshot runner.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXACT COMMAND EXAMPLE:
  python3 scripts/standalone_igv_snapshot_generator.py \\
      --diff-tsv MCBiobank_hg38_comparisons/MCBiobank_hg38_master_report_cleans_common_flags_vs_noflags_diff.tsv \\
      --vironator-dir ./CRAMs \\
      --ref-genome ./REFS/HumanViral_Reference_02-07-2022_modified.renamed.fa \\
      --igv-binary igv.sh \\
      --run-igv
        """
    )
    parser.add_argument(
        "--diff-tsv",
        required=True,
        help="Path to flags vs noflags diff TSV file (sample, virus, flags-noflags, difference)"
    )
    parser.add_argument(
        "--vironator-dir",
        required=True,
        help="Path to directory containing per-sample vironator results"
    )
    parser.add_argument(
        "--ref-genome",
        required=True,
        help="Path to combined multi-FASTA reference genome file"
    )
    parser.add_argument(
        "--output-snapshots-dir",
        default=None,
        help="Optional destination folder for PNG snapshots (defaults to per-sample vironator subfolders)"
    )
    parser.add_argument(
        "--cram-flags",
        default="exogeneSR_viral_clean_filtered.sorted.flags.cram",
        help="Target flags CRAM filename"
    )
    parser.add_argument(
        "--cram-noflags",
        default="exogeneSR_viral_clean_filtered.sorted.cram",
        help="Target noflags CRAM filename"
    )
    parser.add_argument(
        "--cram-additional",
        default="exogeneSR_viral_clean_filtered.sorted.flags.additional.cram",
        help="Target additional difference CRAM filename"
    )
    parser.add_argument(
        "--igv-binary",
        default="igv.sh",
        help="Executable IGV binary command (e.g. igv.sh, igv, xvfb-run -a igv.sh)"
    )
    parser.add_argument(
        "--run-igv",
        action="store_true",
        help="If set, automatically launches IGV batch execution"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    if not os.path.exists(args.diff_tsv):
        print(f"[ERROR] Difference TSV not found at: {args.diff_tsv}")
        sys.exit(1)

    if not os.path.exists(args.ref_genome):
        print(f"[ERROR] Combined reference FASTA not found at: {args.ref_genome}")
        sys.exit(1)

    print(f"[INFO] Reading difference TSV: {args.diff_tsv}")
    df = pd.read_csv(args.diff_tsv, sep="\t")

    # Filter for non-zero difference pairs
    if "difference" in df.columns:
        diff_df = df[df["difference"] != 0].copy()
    elif "flags-noflags" in df.columns:
        diff_df = df[df["flags-noflags"] != 0].copy()
    else:
        diff_df = df.copy()

    print(f"[INFO] Pairs with non-zero differences to visualize: {len(diff_df)}")
    if len(diff_df) == 0:
        print("[INFO] No non-zero difference pairs found. Exiting.")
        sys.exit(0)

    # Batch script output path
    batch_script_path = os.path.abspath("igv_snapshots_batch.igv")
    batch_lines = [
        "new",
        f"genome {os.path.abspath(args.ref_genome)}"
    ]

    valid_pairs_count = 0
    for idx, row in diff_df.iterrows():
        sample = str(row["sample"]).strip()
        virus = str(row["virus"]).strip() if "virus" in row else str(row.get("virus_accession", "")).strip()

        # Candidate directory paths for sample CRAM files
        sample_dir = os.path.abspath(os.path.join(args.vironator_dir, sample))
        flat_dir = os.path.abspath(args.vironator_dir)

        # Check candidate locations for CRAM files (subfolder vs flat)
        tracks_to_load = []
        for label, cram_name in [("flags", args.cram_flags), ("noflags", args.cram_noflags), ("additional", args.cram_additional)]:
            cand_paths = [
                os.path.join(sample_dir, cram_name),
                os.path.join(flat_dir, f"{sample}.{cram_name}"),
                os.path.join(flat_dir, f"{sample}_{cram_name}"),
                os.path.join(flat_dir, cram_name)
            ]
            found = False
            for p in cand_paths:
                if os.path.exists(p) and os.path.getsize(p) > 0:
                    tracks_to_load.append(p)
                    found = True
                    break
            if not found:
                print(f"  [WARNING] Track {label} missing for sample {sample} under {sample_dir} or {flat_dir}")

        if not tracks_to_load:
            print(f"  [SKIP] No CRAM tracks found for sample {sample}")
            continue

        target_snap_dir = os.path.abspath(args.output_snapshots_dir) if args.output_snapshots_dir else sample_dir
        os.makedirs(target_snap_dir, exist_ok=True)

        batch_lines.append("new")
        batch_lines.append(f"snapshotDirectory {target_snap_dir}")
        for track_path in tracks_to_load:
            batch_lines.append(f"load {track_path}")

        batch_lines.append(f"goto {virus}")
        snapshot_filename = f"{sample}_{virus}_igv.png"
        batch_lines.append(f"snapshot {snapshot_filename}")
        valid_pairs_count += 1

    batch_lines.append("exit")

    with open(batch_script_path, "w") as f:
        f.write("\n".join(batch_lines) + "\n")

    print(f"\n[SUCCESS] Saved IGV batch script ({valid_pairs_count} snapshot commands) to: {batch_script_path}")

    if args.run_igv:
        print(f"\n[INFO] Launching IGV with binary: {args.igv_binary}")
        igv_cmd = f"{args.igv_binary} -b \"{batch_script_path}\""
        try:
            res = subprocess.run(igv_cmd, shell=True, check=True)
            print("[SUCCESS] IGV snapshot generation completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] IGV execution failed with return code {e.returncode}.")
            print(f"[HINT] You can run the generated batch script manually in IGV GUI: Tools -> Run Batch Script... -> {batch_script_path}")


if __name__ == "__main__":
    main()
