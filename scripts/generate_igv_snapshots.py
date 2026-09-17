#!/usr/bin/env python3
"""
generate_igv_snapshots.py

Generates IGV batch script files (.igv) and optionally runs IGV in batch mode
to produce PNG snapshots for all sample-virus pairs with non-zero read differences
(difference != 0) in the flags vs noflags comparison TSV.

Visualizes 3 CRAM tracks per sample-virus pair:
1. flags      : exogeneSR_viral_clean_filtered.sorted.flags.cram
2. noflags    : exogeneSR_viral_clean_filtered.sorted.cram
3. additional : exogeneSR_viral_clean_filtered.sorted.flags.additional.cram

Displays the full length of the target virus accession (goto <virus_accession>) without zooming in.
"""

import argparse
import os
import sys
import subprocess
import glob
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate IGV batch scripts and automated PNG snapshots for difference sample-virus pairs."
    )
    parser.add_argument(
        "--diff-tsv",
        default=None,
        help="Path to flags vs noflags diff TSV file (sample, virus, flags-noflags, difference)"
    )
    parser.add_argument(
        "--output-dir",
        default="/mnt/disks/staff",
        help="Base output directory containing vironator results"
    )
    parser.add_argument(
        "--vironator-dirname",
        default="MCBiobank_hg38_vironator",
        help="Vironator output directory name"
    )
    parser.add_argument(
        "--comparison-dirname",
        default="MCBiobank_hg38_comparisons",
        help="Comparison output directory name"
    )
    parser.add_argument(
        "--snapshots-dirname",
        default="MCBiobank_hg38_igv_snapshots",
        help="Directory name to store generated IGV snapshots"
    )
    parser.add_argument(
        "--ref-genome",
        default=None,
        help="Combined multi-FASTA reference genome path"
    )
    parser.add_argument(
        "--dataset",
        default="MCBiobank",
        help="Dataset identifier"
    )
    parser.add_argument(
        "--genome-build",
        default="hg38",
        help="Genome build identifier"
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
        help="Target additional CRAM filename"
    )
    parser.add_argument(
        "--igv-binary",
        default="igv",
        help="Path or binary name for IGV (e.g. igv, igv.sh)"
    )
    parser.add_argument(
        "--run-igv",
        action="store_true",
        help="If set, automatically runs IGV in batch mode after generating batch script"
    )
    parser.add_argument(
        "--phase",
        default="",
        help="Phase identifier (optional)"
    )
    parser.add_argument(
        "--project",
        default="",
        help="Project identifier (optional)"
    )
    parser.add_argument(
        "--output-bucket",
        default=None,
        help="Google Cloud Storage output bucket name (optional)"
    )
    return parser.parse_args()


def resolve_file_path(base_dir, sub_dir, filename):
    candidates = [
        os.path.join(base_dir, sub_dir, filename),
        os.path.join(sub_dir, filename),
        filename
    ]
    for c in candidates:
        if os.path.exists(c):
            return os.path.abspath(c)
    return os.path.abspath(os.path.join(base_dir, sub_dir, filename))


def find_sample_dir(sample, output_dir, vironator_dirname, phase="", project=""):
    """
    Locates the sample vironator directory locally or under output_dir ({dataset}_{build}_vironator/{sample}).
    """
    candidates = [
        os.path.join(output_dir, vironator_dirname, sample),
        os.path.join(".", vironator_dirname, sample),
        os.path.join(output_dir, sample),
        os.path.join(".", sample)
    ]
    for c in candidates:
        if os.path.isdir(c):
            return os.path.abspath(c)
    return os.path.abspath(os.path.join(".", vironator_dirname, sample))


def fetch_cram_from_gcs(bucket, vironator_dirname, sample, cram_filename, target_dir, phase="", project=""):
    """
    Pulls missing CRAM file from GCS bucket (gs://{bucket}/{vironator_dirname}/{sample}/{cram_filename}).
    """
    os.makedirs(target_dir, exist_ok=True)
    target_path = os.path.join(target_dir, cram_filename)
    if os.path.exists(target_path) and os.path.getsize(target_path) > 0:
        return target_path

    gcs_candidates = [
        f"gs://{bucket}/{vironator_dirname}/{sample}/{cram_filename}",
        f"gs://{bucket}/{sample}/{cram_filename}"
    ]

    for gcs_src in gcs_candidates:
        gcs_crai = f"{gcs_src}.crai"
        check_cmd = f"gsutil -q stat \"{gcs_src}\" 2>/dev/null"
        if subprocess.run(check_cmd, shell=True).returncode == 0:
            print(f"  [GCS FETCH] Pulling CRAM from GCS: {gcs_src} -> {target_path}")
            cmd = f"gsutil -q cp \"{gcs_src}\" \"{target_path}\" 2>/dev/null && gsutil -q cp \"{gcs_crai}\" \"{target_path}.crai\" 2>/dev/null || true"
            subprocess.run(cmd, shell=True)

            if os.path.exists(target_path) and os.path.getsize(target_path) > 0:
                return target_path

    return None



def main():
    args = parse_args()
    prefix = f"{args.dataset}_{args.genome_build}"

    # Determine diff TSV location
    diff_tsv_path = args.diff_tsv
    if not diff_tsv_path:
        default_diff_name = f"{prefix}_master_report_cleans_common_flags_vs_noflags_diff.tsv"
        diff_tsv_path = resolve_file_path(args.output_dir, args.comparison_dirname, default_diff_name)

    print(f"[INFO] Using difference TSV: {diff_tsv_path}")
    if not os.path.exists(diff_tsv_path):
        print(f"[ERROR] Difference TSV not found at: {diff_tsv_path}")
        sys.exit(1)

    # Read diff TSV
    df = pd.read_csv(diff_tsv_path, sep="\t")
    print(f"[INFO] Total pairs in difference TSV: {len(df)}")

    # Filter for difference != 0
    # Columns expected: sample, virus (or virus_accession), flags-noflags (or difference_str), difference
    if "difference" in df.columns:
        diff_df = df[df["difference"] != 0].copy()
    elif "flags-noflags" in df.columns:
        diff_df = df[df["flags-noflags"] != 0].copy()
    else:
        # Fallback: take all rows if difference column missing
        diff_df = df.copy()

    print(f"[INFO] Pairs with non-zero differences: {len(diff_df)}")

    if len(diff_df) == 0:
        print("[INFO] No non-zero difference pairs to visualize. Exiting.")
        sys.exit(0)

    # Prepare output snapshot directory (fallback to local work_dir when mounted disk is unavailable)
    snapshots_dir = os.path.join(args.output_dir, args.snapshots_dirname)
    try:
        os.makedirs(snapshots_dir, exist_ok=True)
    except Exception:
        snapshots_dir = os.path.abspath(os.path.join(".", args.snapshots_dirname))
        os.makedirs(snapshots_dir, exist_ok=True)
    print(f"[INFO] Output snapshots directory: {snapshots_dir}")

    # Resolve reference genome
    ref_genome = args.ref_genome
    if not ref_genome:
        ref_genome = "/mnt/disks/staff/refs/HumanViral_Reference_02-07-2022_SnapGene_plasmids_modified_mm39_ms_modified.fa"

    # Prepare IGV batch script file
    batch_script_path = os.path.join(snapshots_dir, f"{prefix}_igv_snapshots.igv")
    
    batch_lines = [
        "new",
        f"genome {ref_genome}"
    ]

    valid_pairs_count = 0
    for idx, row in diff_df.iterrows():
        sample = str(row["sample"]).strip()
        virus = str(row["virus"]).strip() if "virus" in row else str(row.get("virus_accession", "")).strip()

        sample_dir = find_sample_dir(sample, args.output_dir, args.vironator_dirname, phase=args.phase, project=args.project)
        
        cram_flags_path = os.path.join(sample_dir, args.cram_flags)
        cram_noflags_path = os.path.join(sample_dir, args.cram_noflags)
        cram_add_path = os.path.join(sample_dir, args.cram_additional)

        # Check track availability and fetch from GCS bucket if missing locally
        tracks_to_load = []
        for label, filename, path in [
            ("flags", args.cram_flags, cram_flags_path),
            ("noflags", args.cram_noflags, cram_noflags_path),
            ("additional", args.cram_additional, cram_add_path)
        ]:
            if os.path.exists(path) and os.path.getsize(path) > 0:
                tracks_to_load.append(path)
            elif args.output_bucket:
                fetched_path = fetch_cram_from_gcs(args.output_bucket, args.vironator_dirname, sample, filename, sample_dir, phase=args.phase, project=args.project)
                if fetched_path:
                    tracks_to_load.append(fetched_path)
                else:
                    print(f"  [WARNING] Track {label} missing locally and on GCS for sample {sample}: {filename}")
            else:
                print(f"  [WARNING] Track {label} missing for sample {sample}: {path}")

        if not tracks_to_load:
            print(f"  [SKIP] No CRAM tracks found for sample {sample} under {sample_dir}")
            continue

        batch_lines.append("new")
        batch_lines.append(f"snapshotDirectory {sample_dir}")
        for track_path in tracks_to_load:
            batch_lines.append(f"load {track_path}")

        # Navigate to virus accession full length
        batch_lines.append(f"goto {virus}")
        
        # Save snapshot inside sample's vironator directory
        snapshot_filename = f"{sample}_{virus}_igv.png"
        batch_lines.append(f"snapshot {snapshot_filename}")
        valid_pairs_count += 1

    batch_lines.append("exit")

    # Write IGV batch script
    with open(batch_script_path, "w") as f:
        f.write("\n".join(batch_lines) + "\n")

    print(f"[SUCCESS] Wrote IGV batch script ({valid_pairs_count} snapshot commands) to: {batch_script_path}")

    # Optionally execute IGV batch mode or igv-reports
    if args.run_igv:
        igv_bin = args.igv_binary
        if igv_bin.startswith("xvfb-run"):
            check_xvfb = subprocess.run("which xvfb-run 2>/dev/null", shell=True).returncode
            if check_xvfb != 0:
                print("  [WARNING] 'xvfb-run' not found on PATH. Falling back to direct IGV invocation ('igv.sh').")
                igv_bin = igv_bin.split()[-1]

        # Check if igv-reports is requested or available
        if "igv-reports" in igv_bin or "create_report" in igv_bin:
            print(f"[INFO] Launching IGV-Reports HTML generator for {valid_pairs_count} pairs...")
            for idx, row in diff_df.iterrows():
                sample = str(row["sample"]).strip()
                virus = str(row["virus"]).strip() if "virus" in row else str(row.get("virus_accession", "")).strip()
                sample_dir = find_sample_dir(sample, args.output_dir, args.vironator_dirname, phase=args.phase, project=args.project)
                cram_flags_path = os.path.join(sample_dir, args.cram_flags)
                cram_noflags_path = os.path.join(sample_dir, args.cram_noflags)
                cram_add_path = os.path.join(sample_dir, args.cram_additional)
                tracks = [p for p in [cram_flags_path, cram_noflags_path, cram_add_path] if os.path.exists(p)]
                if tracks:
                    html_out = os.path.join(sample_dir, f"{sample}_{virus}_igv_report.html")
                    cmd = f"create_report {virus} --fasta \"{ref_genome}\" --tracks {' '.join(tracks)} --output \"{html_out}\" 2>/dev/null || true"
                    subprocess.run(cmd, shell=True)
            print("[SUCCESS] Completed IGV HTML report generation.")
        else:
            print(f"[INFO] Launching IGV batch mode with binary: {igv_bin}")
            igv_cmd = f"{igv_bin} -b {batch_script_path}"
            try:
                res = subprocess.run(igv_cmd, shell=True, check=True, capture_output=True, text=True)
                print("[INFO] IGV stdout:\n", res.stdout)
            except subprocess.CalledProcessError as e:
                print(f"[WARNING] IGV batch GUI execution failed (likely headless X11 missing): {e.stderr}")
                print("[NOTE] Batch script is preserved at: ", batch_script_path)
                print("[NOTE] You can run 'create_report' (igv-reports) or load the .igv script directly in desktop IGV GUI.")


if __name__ == "__main__":
    main()
