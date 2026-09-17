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
        "--ref-genome",
        default=None,
        help="Reference FASTA path for CRAM decoding/encoding (optional)"
    )
    return parser.parse_args()


def find_sample_vironator_dir(sample_id, base_dir, vironator_dirname):
    """
    Locates the specific sample directory inside the vironator output structure.
    """
    vironator_base = os.path.join(base_dir, vironator_dirname)
    if not os.path.exists(vironator_base):
        vironator_base = base_dir

    # Search pattern for sample directory
    pattern = os.path.join(vironator_base, "**", sample_id)
    matches = [d for d in glob.glob(pattern, recursive=True) if os.path.isdir(d)]
    
    if matches:
        return matches[0]
    
    # Fallback flat path
    flat_path = os.path.join(vironator_base, sample_id)
    return flat_path if os.path.exists(flat_path) else None


def process_master_report_comparison(master_report_path, output_stats_dir):
    """
    Implements the complete workflow described in OLD/MCBiobank_flags_vs_noflags_workflow.docx:
    1. Filter out 'None' rows.
    2. Split into noflags vs flags entries (based on Source_File / CRAM strategy or filename).
    3. Compare sample-virus pairs (Sample_ID + Virus_Accession) to generate:
       - common list (shared sample-virus pairs with read count comparison)
       - noflags unique list
       - flags unique list (and summary of extra hits per virus)
    4. Generate summary TSVs / mini-reports in output_stats_dir.
    5. Return sample list for difference CRAM generation.
    """
    if not os.path.exists(master_report_path):
        local_repo_master = os.path.basename(master_report_path)
        if os.path.exists(local_repo_master):
            master_report_path = local_repo_master
        else:
            print(f"[ERROR] Master report TSV not found at '{master_report_path}'.")
            print("[IMPORTANT] Please ensure the master report TSV is copied directly to your cloned VIROnator repository directory.")
            sys.exit(1)

    print(f"[INFO] Reading master report: {master_report_path}")
    df = pd.read_csv(master_report_path, sep="\t", dtype=str)
    df.columns = [c.strip() for c in df.columns]

    # Filter out rows containing 'None' or empty values
    df_clean = df[~df.isin(["None", "none", None]).any(axis=1)].copy()

    # Identify key columns
    sample_col = df.columns[0]  # Sample_ID
    virus_acc_col = df.columns[1] if len(df.columns) > 1 else df.columns[0]  # Virus_Accession
    virus_name_col = df.columns[9] if len(df.columns) > 9 else virus_acc_col
    reads_col = df.columns[3] if len(df.columns) > 3 else df.columns[0]
    source_col = [c for c in df.columns if "source" in c.lower() or "file" in c.lower()]
    source_col = source_col[0] if source_col else df.columns[-1]

    # Split into flags and noflags DataFrames
    is_flags = df_clean[source_col].str.contains("flags", case=False, na=False)
    df_flags = df_clean[is_flags].copy()
    df_noflags = df_clean[~is_flags].copy()

    # Ensure output stats dir exists
    os.makedirs(output_stats_dir, exist_ok=True)

    # Save cleaned splits
    noflags_tsv = os.path.join(output_stats_dir, "master_report_cleans_noflags.tsv")
    flags_tsv = os.path.join(output_stats_dir, "master_report_cleans_flags.tsv")
    df_noflags.to_csv(noflags_tsv, sep="\t", index=False)
    df_flags.to_csv(flags_tsv, sep="\t", index=False)

    # Build comparison keys (Sample_ID + Virus_Accession)
    df_flags['key'] = df_flags[sample_col] + "___" + df_flags[virus_acc_col]
    df_noflags['key'] = df_noflags[sample_col] + "___" + df_noflags[virus_acc_col]

    flags_keys = set(df_flags['key'])
    noflags_keys = set(df_noflags['key'])

    common_keys = flags_keys.intersection(noflags_keys)
    flags_unique_keys = flags_keys - noflags_keys
    noflags_unique_keys = noflags_keys - flags_keys

    # 1. Output Unique Lists (Mini-Reports)
    df_flags_unique = df_flags[df_flags['key'].isin(flags_unique_keys)].drop(columns=['key'])
    df_noflags_unique = df_noflags[df_noflags['key'].isin(noflags_unique_keys)].drop(columns=['key'])

    flags_unique_path = os.path.join(output_stats_dir, "master_report_cleans_flags_unique_list.tsv")
    noflags_unique_path = os.path.join(output_stats_dir, "master_report_cleans_noflags_unique_list.tsv")
    df_flags_unique.to_csv(flags_unique_path, sep="\t", index=False)
    df_noflags_unique.to_csv(noflags_unique_path, sep="\t", index=False)

    # 2. Summarize Flags-Unique Hits by Virus
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
        summary_path = os.path.join(output_stats_dir, "master_report_cleans_flags_unique_list_summary.tsv")
        df_summary.to_csv(summary_path, sep="\t", index=False)
        print(f"[REPORT] Saved flags-unique viral summary: {summary_path}")

    # 3. Read Count Differences for Shared Hits
    diff_records = []
    flags_map = df_flags.set_index('key')
    noflags_map = df_noflags.set_index('key')

    for key in sorted(common_keys):
        r_flag = flags_map.loc[key]
        r_noflag = noflags_map.loc[key]

        # Handle duplicates if any by taking first or sum
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
    diff_path = os.path.join(output_stats_dir, "master_report_cleans_common_flags_vs_noflags_diff.tsv")
    df_diff.to_csv(diff_path, sep="\t", index=False)
    print(f"[REPORT] Saved common hits read count difference TSV: {diff_path}")

    samples = sorted(df[sample_col].dropna().str.strip().unique())
    return samples, df_diff


def main():
    args = parse_args()

    print("======================================================================")
    print("VIROnator Flag Difference CRAM Generator & Mini-Report Module")
    print("======================================================================")

    samples = []
    if args.input_tsv:
        samples, _ = parse_custom_input_tsv(args.input_tsv)
    elif args.master_report:
        stats_dir = os.path.join(args.output_dir, "MCBiobank_hg38_stats")
        samples, _ = process_master_report_comparison(args.master_report, stats_dir)
    else:
        repo_master = "MCBiobank_hg38_master_report.tsv"
        if os.path.exists(repo_master):
            stats_dir = os.path.join(args.output_dir, "MCBiobank_hg38_stats")
            samples, _ = process_master_report_comparison(repo_master, stats_dir)
        else:
            print("[ERROR] Neither --input-tsv nor --master-report was provided.")
            print("[NOTE] Please provide an input TSV file or ensure the master report TSV is placed in the VIROnator repository directory.")
            sys.exit(1)

    total_samples = len(samples)
    print(f"\n[START] Processing CRAM differences for {total_samples} samples sequentially...\n")

    success_count = 0
    for idx, sample_id in enumerate(samples, start=1):
        print(f"[{idx}/{total_samples}] Processing sample: {sample_id}")
        sample_dir = find_sample_vironator_dir(sample_id, args.output_dir, args.vironator_dirname)

        if not sample_dir or not os.path.exists(sample_dir):
            print(f"  [WARNING] Vironator directory not found for sample {sample_id}. Skipping.")
            continue

        ok = generate_difference_cram_for_sample(sample_id, sample_dir, args)
        if ok:
            success_count += 1

    print("\n======================================================================")
    print(f"[COMPLETED] Successfully processed {success_count}/{total_samples} samples.")
    print("======================================================================")


if __name__ == "__main__":
    main()

