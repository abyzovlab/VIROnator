#!/usr/bin/env python3
"""
generate_flag_difference_crams.py

Generates exogeneSR_viral_clean_filtered.sorted.flags.additional.cram for each sample by:
1. Checking if a custom difference input TSV (sample, virus, flags-noflags, difference) is provided.
   - If NOT provided, extracts samples/differences from the master report TSV (cloned in the local VIROnator repo directory).
2. Iterating sequentially over samples to locate their vironator directories and the two target CRAM files:
   - File 1 (noflags): exogeneSR_viral_clean_filtered.sorted.cram
   - File 2 (flags)  : exogeneSR_viral_clean_filtered.sorted.flags.cram
3. Subtracting noflags reads from flags reads to create:
   - exogeneSR_viral_clean_filtered.sorted.flags.additional.cram
4. Outputting real-time progress to stdout.
"""

import argparse
import os
import sys
import subprocess
import glob
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate flag difference CRAM files per sample without batch jobs."
    )
    parser.add_argument(
        "--input-tsv",
        default=None,
        help="Path to custom 4-column difference TSV file (sample, virus, flags-noflags, difference)"
    )
    parser.add_argument(
        "--master-report",
        default=None,
        help="Path to master report TSV in VIROnator repo directory (used if input-tsv is not provided)"
    )
    parser.add_argument(
        "--cram-noflags",
        default="exogeneSR_viral_clean_filtered.sorted.cram",
        help="Target noflags CRAM filename"
    )
    parser.add_argument(
        "--cram-flags",
        default="exogeneSR_viral_clean_filtered.sorted.flags.cram",
        help="Target flags CRAM filename"
    )
    parser.add_argument(
        "--cram-additional",
        default="exogeneSR_viral_clean_filtered.sorted.flags.additional.cram",
        help="Output additional difference CRAM filename"
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
        "--stats-dirname",
        default="MCBiobank_hg38_stats",
        help="Stats output directory name"
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
    parser.add_argument(
        "--ref-genome",
        default=None,
        help="Reference FASTA path for CRAM decoding/encoding (optional)"
    )
    parser.add_argument(
        "--generate-igv-snapshots",
        action="store_true",
        help="If set, generates IGV batch script and PNG snapshots for non-zero difference pairs"
    )
    parser.add_argument(
        "--igv-binary",
        default="igv",
        help="Path or binary name for IGV (e.g. igv, igv.sh)"
    )
    parser.add_argument(
        "--ref-viral-fasta",
        default=None,
        help="Combined multi-FASTA viral reference genome path for IGV"
    )
    return parser.parse_args()


def find_sample_vironator_dir(sample_id, base_dir, vironator_dirname, gcs_bucket=None, phase="", project=""):
    """
    Locates and fetches sample CRAM directory directly from GCS output bucket using output_bucket, dataset, build, phase, project variables.
    """
    local_sample_dir = os.path.join(".", vironator_dirname, sample_id)

    # If already downloaded locally in previous step, reuse
    if os.path.exists(local_sample_dir) and os.path.isdir(local_sample_dir):
        return local_sample_dir

    if gcs_bucket:
        # Candidate GCS paths (flat, phase/project structured)
        phase_part = f"phase{phase}/" if phase and str(phase).strip() and str(phase).strip().lower() not in ["none", "0"] else ""
        project_part = f"{project}/" if project and str(project).strip() and str(project).strip().lower() not in ["none", "0", "base"] else ""

        gcs_candidates = [
            f"gs://{gcs_bucket}/{vironator_dirname}/{phase_part}{project_part}{sample_id}/",
            f"gs://{gcs_bucket}/{vironator_dirname}/{sample_id}/",
            f"gs://{gcs_bucket}/{sample_id}/"
        ]

        for gcs_sample_path in gcs_candidates:
            check_cmd = f"gsutil -q stat \"{gcs_sample_path}*\" 2>/dev/null"
            if subprocess.run(check_cmd, shell=True).returncode == 0:
                print(f"  [GCS FETCH] Pulling CRAM files from GCS: {gcs_sample_path} -> {local_sample_dir}", flush=True)
                os.makedirs(local_sample_dir, exist_ok=True)
                dl_cmd = f"gsutil -q -m cp -r \"{gcs_sample_path}*\" \"{local_sample_dir}/\" 2>/dev/null || true"
                subprocess.run(dl_cmd, shell=True)

                if os.path.exists(local_sample_dir) and len(os.listdir(local_sample_dir)) > 0:
                    return local_sample_dir

    # Fallback check for any existing local path
    for cand in [os.path.join(base_dir, vironator_dirname, sample_id), os.path.join(".", sample_id)]:
        if os.path.exists(cand) and os.path.isdir(cand):
            return cand

    return None


def parse_custom_input_tsv(input_tsv_path):
    """
    Parses custom 4-column difference TSV (sample, virus, flags-noflags, difference).
    """
    if not os.path.exists(input_tsv_path):
        print(f"[ERROR] Custom input TSV file not found: {input_tsv_path}", flush=True)
        sys.exit(1)

    print(f"[INFO] Reading custom difference TSV: {input_tsv_path}", flush=True)
    df = pd.read_csv(input_tsv_path, sep="\t", dtype=str)
    sample_col = df.columns[0]
    samples = sorted(df[sample_col].dropna().str.strip().unique())
    return samples, df


def generate_difference_cram_for_sample(sample_id, sample_dir, args):
    """
    Generates additional CRAM file for a single sample by finding read IDs present in flags CRAM
    that are absent from noflags CRAM.
    """
    cram_noflags_path = os.path.join(sample_dir, args.cram_noflags)
    cram_flags_path = os.path.join(sample_dir, args.cram_flags)
    cram_add_path = os.path.join(sample_dir, args.cram_additional)

    if not os.path.exists(cram_flags_path):
        print(f"  [SKIP] Flags CRAM missing for sample {sample_id}: {cram_flags_path}", flush=True)
        return False

    if not os.path.exists(cram_noflags_path):
        print(f"  [SKIP] Noflags CRAM missing for sample {sample_id}: {cram_noflags_path}", flush=True)
        return False

    ref_flag = f"-T \"{args.ref_genome}\"" if args.ref_genome and os.path.exists(args.ref_genome) else ""

    # Shell pipeline to extract difference read names and filter flags CRAM
    tmp_noflags_ids = os.path.join(sample_dir, ".tmp_noflags_ids.txt")
    tmp_flags_ids = os.path.join(sample_dir, ".tmp_flags_ids.txt")
    tmp_diff_ids = os.path.join(sample_dir, ".tmp_diff_ids.txt")

    try:
        # Extract sorted QNAMEs
        cmd_noflags = f"samtools view {ref_flag} \"{cram_noflags_path}\" | cut -f1 | sort -u > \"{tmp_noflags_ids}\""
        cmd_flags = f"samtools view {ref_flag} \"{cram_flags_path}\" | cut -f1 | sort -u > \"{tmp_flags_ids}\""
        subprocess.run(cmd_noflags, shell=True, check=True)
        subprocess.run(cmd_flags, shell=True, check=True)

        # Comm difference
        cmd_comm = f"comm -13 \"{tmp_noflags_ids}\" \"{tmp_flags_ids}\" > \"{tmp_diff_ids}\""
        subprocess.run(cmd_comm, shell=True, check=True)

        diff_count = 0
        if os.path.exists(tmp_diff_ids):
            with open(tmp_diff_ids, "r") as f:
                diff_count = sum(1 for _ in f)

        if diff_count == 0:
            print(f"  [INFO] 0 difference reads for sample {sample_id}. Creating empty header CRAM.", flush=True)
            cmd_empty = f"samtools view -H {ref_flag} \"{cram_flags_path}\" | samtools view -b - | samtools convert {ref_flag} -O cram -o \"{cram_add_path}\""
            subprocess.run(cmd_empty, shell=True, check=True)
        else:
            print(f"  [BUILD] Extracting {diff_count} additional reads into {args.cram_additional}...", flush=True)
            cmd_filter = f"samtools view -N \"{tmp_diff_ids}\" {ref_flag} -O cram -o \"{cram_add_path}\" \"{cram_flags_path}\""
            subprocess.run(cmd_filter, shell=True, check=True)

        # Clean temp files
        for tmp_f in [tmp_noflags_ids, tmp_flags_ids, tmp_diff_ids]:
            if os.path.exists(tmp_f):
                os.remove(tmp_f)

        print(f"  [SUCCESS] Created: {cram_add_path}", flush=True)

        # Upload generated additional CRAM back to GCS bucket if output_bucket is specified
        if args.output_bucket:
            phase_part = f"phase{args.phase}/" if args.phase and str(args.phase).strip() and str(args.phase).strip().lower() not in ["none", "0"] else ""
            project_part = f"{args.project}/" if args.project and str(args.project).strip() and str(args.project).strip().lower() not in ["none", "0", "base"] else ""
            gcs_dest_sample = f"gs://{args.output_bucket}/{args.vironator_dirname}/{phase_part}{project_part}{sample_id}/"
            
            print(f"  [GCS UPLOAD] Syncing generated difference CRAM to {gcs_dest_sample}", flush=True)
            upload_cmd = f"gsutil -q cp \"{cram_add_path}\" \"{gcs_dest_sample}\" 2>/dev/null || true"
            subprocess.run(upload_cmd, shell=True)

        # Cleanup local staging sample folder if it was staged under local ./vironator_dirname
        local_staging_dir = os.path.abspath(os.path.join(".", args.vironator_dirname, sample_id))
        if os.path.abspath(sample_dir) == local_staging_dir and os.path.exists(local_staging_dir):
            import shutil
            shutil.rmtree(local_staging_dir, ignore_errors=True)
            print(f"  [CLEANUP] Removed temporary local sample staging directory: {local_staging_dir}", flush=True)

        return True

    except Exception as e:
        print(f"  [ERROR] Failed generating difference CRAM for {sample_id}: {e}", flush=True)
        for tmp_f in [tmp_noflags_ids, tmp_flags_ids, tmp_diff_ids]:
            if os.path.exists(tmp_f):
                try: os.remove(tmp_f)
                except Exception: pass
        return False


def process_master_report_comparison(master_report_path, output_stats_dir, dataset="MCBiobank", genome_build="hg38", phase="", project=""):
    """
    Implements the complete workflow described in OLD/MCBiobank_flags_vs_noflags_workflow.docx:
    Outputs TSV files prepended with prefix: {dataset}_{genome_build}[_{phase}][_{project}]
    """
    if not os.path.exists(master_report_path):
        local_repo_master = os.path.basename(master_report_path)
        if os.path.exists(local_repo_master):
            master_report_path = local_repo_master
        else:
            print(f"[ERROR] Master report TSV not found at '{master_report_path}'.")
            print("[IMPORTANT] Please ensure the master report TSV is copied directly to your cloned VIROnator repository directory.")
            sys.exit(1)

    # Build dynamic filename prefix
    parts = [str(dataset).strip(), str(genome_build).strip()]
    if phase and str(phase).strip() and str(phase).strip().lower() not in ["none", "0"]:
        parts.append(str(phase).strip())
    if project and str(project).strip() and str(project).strip().lower() not in ["none", "0", "base"]:
        parts.append(str(project).strip())
    prefix = "_".join(parts)

    print(f"[INFO] Reading master report: {master_report_path}")
    df = pd.read_csv(master_report_path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # Identify key columns
    sample_col = df.columns[0]  # Sample_ID
    virus_acc_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]  # Virus_Accession
    virus_name_col = df.columns[9] if len(df.columns) > 9 else virus_acc_col
    reads_col = df.columns[3] if len(df.columns) > 3 else df.columns[0]
    source_col = [c for c in df.columns if "source" in c.lower() or "file" in c.lower()]
    source_col = source_col[0] if source_col else df.columns[-1]

    # Step 0: Filter out rows containing 'None', empty values, or <= 0 mapped reads
    df_clean = df[~df.isin(["None", "none", "NONE", None]).any(axis=1)].copy()
    if reads_col in df_clean.columns:
        df_clean[reads_col] = pd.to_numeric(df_clean[reads_col], errors='coerce').fillna(0)
        df_clean = df_clean[df_clean[reads_col] > 0].copy()

    # Step 0 Split: flags vs noflags matching target CRAM filename strings
    flags_target = "exogeneSR_viral_clean_filtered.sorted.flags"
    noflags_target = "exogeneSR_viral_clean_filtered.sorted.cram"

    is_flags = df_clean[source_col].str.contains(flags_target, case=False, na=False)
    # noflags matches exogeneSR_viral_clean_filtered.sorted.cram but excludes flags
    is_noflags = df_clean[source_col].str.contains(noflags_target, case=False, na=False) & (~is_flags)

    df_flags = df_clean[is_flags].copy()
    df_noflags = df_clean[is_noflags].copy()

    # Ensure output stats dir exists
    os.makedirs(output_stats_dir, exist_ok=True)

    # Save cleaned splits per Step 0 with dynamic prefix
    cleans_tsv = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans.tsv")
    noflags_tsv = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_noflags.tsv")
    flags_tsv = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_flags.tsv")

    df_clean.to_csv(cleans_tsv, sep="\t", index=False)
    df_noflags.to_csv(noflags_tsv, sep="\t", index=False)
    df_flags.to_csv(flags_tsv, sep="\t", index=False)

    # Build comparison keys (Sample_ID + Virus_Accession)
    df_flags['key'] = df_flags[sample_col].astype(str) + "___" + df_flags[virus_acc_col].astype(str)
    df_noflags['key'] = df_noflags[sample_col].astype(str) + "___" + df_noflags[virus_acc_col].astype(str)

    flags_keys = set(df_flags['key'].dropna())
    noflags_keys = set(df_noflags['key'].dropna())

    common_keys = flags_keys.intersection(noflags_keys)
    flags_unique_keys = flags_keys - noflags_keys
    noflags_unique_keys = noflags_keys - flags_keys

    # 1. Output Common List (Step 1 & Step 3) and Unique Lists (Step 1)
    df_flags_common = df_flags[df_flags['key'].isin(common_keys)].copy()
    df_noflags_common = df_noflags[df_noflags['key'].isin(common_keys)].copy()

    df_flags_common.insert(0, 'Source_Strategy', 'flags')
    df_noflags_common.insert(0, 'Source_Strategy', 'noflags')

    df_common_concat = pd.concat([df_flags_common, df_noflags_common], ignore_index=True)
    df_common_ordered = df_common_concat.sort_values(by=[sample_col, virus_acc_col, 'Source_Strategy']).drop(columns=['key'])

    common_ordered_path = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_common_list_ordered.tsv")
    df_common_ordered.to_csv(common_ordered_path, sep="\t", index=False)

    df_flags_unique = df_flags[df_flags['key'].isin(flags_unique_keys)].drop(columns=['key'])
    df_noflags_unique = df_noflags[df_noflags['key'].isin(noflags_unique_keys)].drop(columns=['key'])

    flags_unique_path = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_flags_unique_list.tsv")
    noflags_unique_path = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_noflags_unique_list.tsv")
    df_flags_unique.to_csv(flags_unique_path, sep="\t", index=False)
    df_noflags_unique.to_csv(noflags_unique_path, sep="\t", index=False)

    # 2. Summarize Flags-Unique Hits by Virus (Step 2)
    if not df_flags_unique.empty:
        df_flags_unique[reads_col] = pd.to_numeric(df_flags_unique[reads_col], errors='coerce').fillna(0)
        summary_rows = []
        grouped = df_flags_unique.groupby([virus_acc_col, virus_name_col])
        total_unique_rows = len(df_flags_unique)

        for (acc, name), group in grouped:
            num_rows = len(group)
            pct_rows = (num_rows / total_unique_rows) * 100.0 if total_unique_rows > 0 else 0.0
            tot_reads = group[reads_col].sum()
            max_reads = group[reads_col].max()
            summary_rows.append({
                "accession": acc,
                "virus_name": name,
                "number_of_rows": num_rows,
                "percentage_of_rows": round(pct_rows, 2),
                "total_reads": int(tot_reads),
                "max_read_count_observed": int(max_reads)
            })

        df_summary = pd.DataFrame(summary_rows)
        if not df_summary.empty:
            df_summary = df_summary.sort_values(by="total_reads", ascending=False)
        summary_path = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_flags_unique_list_summary.tsv")
        df_summary.to_csv(summary_path, sep="\t", index=False)
        print(f"[REPORT] Saved flags-unique viral summary: {summary_path}")

    # 3. Read Count Differences for Shared Hits (Step 4)
    diff_records = []
    flags_map = df_flags.set_index('key')
    noflags_map = df_noflags.set_index('key')

    for key in sorted(common_keys, key=lambda k: str(k)):
        r_flag = flags_map.loc[key]
        r_noflag = noflags_map.loc[key]

        s_id = r_flag[sample_col].iloc[0] if isinstance(r_flag, pd.DataFrame) else r_flag[sample_col]
        v_acc = r_flag[virus_acc_col].iloc[0] if isinstance(r_flag, pd.DataFrame) else r_flag[virus_acc_col]

        f_reads = pd.to_numeric(r_flag[reads_col], errors='coerce').sum() if isinstance(r_flag, pd.DataFrame) else float(r_flag[reads_col])
        nf_reads = pd.to_numeric(r_noflag[reads_col], errors='coerce').sum() if isinstance(r_noflag, pd.DataFrame) else float(r_noflag[reads_col])
        diff = f_reads - nf_reads

        diff_records.append({
            "sample": s_id,
            "virus": v_acc,
            "flags-noflags": f"{int(f_reads)}-{int(nf_reads)}",
            "difference": int(diff)
        })

    df_diff = pd.DataFrame(diff_records)
    diff_path = os.path.join(output_stats_dir, f"{prefix}_master_report_cleans_common_flags_vs_noflags_diff.tsv")
    df_diff.to_csv(diff_path, sep="\t", index=False)
    print(f"[REPORT] Saved common hits read count difference TSV: {diff_path}")

    # Filter samples to ONLY those with non-zero differences (difference != 0) in common_flags_vs_noflags_diff
    if not df_diff.empty:
        non_zero_diff_df = df_diff[df_diff['difference'] != 0]
        samples = sorted(non_zero_diff_df['sample'].dropna().str.strip().unique())
    else:
        samples = []

    return samples, df_diff


def main():
    args = parse_args()

    print("======================================================================", flush=True)
    print("VIROnator Flag Difference CRAM Generator & Mini-Report Module", flush=True)
    print("======================================================================", flush=True)

    # Resolve stats_dir with fallback to local work_dir when mounted disk is unavailable
    stats_dir = os.path.join(args.output_dir, args.stats_dirname)
    try:
        os.makedirs(stats_dir, exist_ok=True)
    except Exception:
        stats_dir = os.path.join(".", args.stats_dirname)
        os.makedirs(stats_dir, exist_ok=True)

    samples = []
    if args.input_tsv:
        samples, _ = parse_custom_input_tsv(args.input_tsv)
    elif args.master_report:
        samples, _ = process_master_report_comparison(args.master_report, stats_dir, args.dataset, args.genome_build, args.phase, args.project)
    else:
        repo_master = f"{args.dataset}_{args.genome_build}_master_report.tsv"
        if not os.path.exists(repo_master):
            repo_master = "MCBiobank_hg38_master_report.tsv"
        if os.path.exists(repo_master):
            samples, _ = process_master_report_comparison(repo_master, stats_dir, args.dataset, args.genome_build, args.phase, args.project)
        else:
            print("[ERROR] Neither --input-tsv nor --master-report was provided.", flush=True)
            print("[NOTE] Please provide an input TSV file or ensure the master report TSV is placed in the VIROnator repository directory.", flush=True)
            sys.exit(1)

    total_samples = len(samples)
    print(f"\n[START] Processing CRAM differences for {total_samples} samples sequentially...\n", flush=True)

    success_count = 0
    for idx, sample_id in enumerate(samples, start=1):
        print(f"[{idx}/{total_samples}] Processing sample: {sample_id}", flush=True)
        sample_dir = find_sample_vironator_dir(sample_id, args.output_dir, args.vironator_dirname, gcs_bucket=args.output_bucket, phase=args.phase, project=args.project)

        if not sample_dir or not os.path.exists(sample_dir):
            print(f"  [WARNING] Sample directory for '{sample_id}' not found in GCS bucket 'gs://{args.output_bucket}/{args.vironator_dirname}/'. Skipping.", flush=True)
            continue

        ok = generate_difference_cram_for_sample(sample_id, sample_dir, args)
        if ok:
            success_count += 1

    # Cleanup top-level local vironator directory if empty or remaining
    local_vir_parent = os.path.abspath(os.path.join(".", args.vironator_dirname))
    if os.path.exists(local_vir_parent):
        import shutil
        shutil.rmtree(local_vir_parent, ignore_errors=True)
        print(f"[CLEANUP] Removed temporary top-level staging directory: {local_vir_parent}", flush=True)

    print("\n======================================================================", flush=True)
    print(f"[COMPLETED] Successfully processed CRAM differences for {success_count}/{total_samples} samples.", flush=True)
    print("======================================================================", flush=True)


if __name__ == "__main__":
    main()

