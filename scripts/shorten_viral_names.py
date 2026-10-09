#!/usr/bin/env python3
"""
shorten_viral_names.py
Standalone Utility Script for VIROnator

Takes an input viral reference rename map file (.renamed_map.tsv) and automatically
determines the output path in the same directory, appending '_short.tsv'.

Output File Structure (3 Columns):
  Column 1: final_clean_id (Accession ID)
  Column 2: original_header (Unmodified original header)
  Column 3: short_name (Shortened viral name without leading accession ID)

Example Output:
  HE974372 \t HE974372 HBV genotype D4 complete genome, isolate Mart-B36 \t HBV genotype D4 isolate Mart-B36

Rules Enforced:
1. Substring Removals:
   - " genomic RNA, complete genome, strain:"
   - ", complete genome, strain:"
   - " genomic RNA, complete"
   - " genomic DNA, partial"
   - " genomic DNA."
   - ", complete genome."
   - ", complete genome"
   - "sense strand"

2. Abbreviation Replacements:
   - "Torque teno mini virus" -> "TTV mini"
   - "Torque teno virus"      -> "TTV"
   - "Simian adenovirus"      -> "Simian Ad"
   - "Human papillomavirus"   -> "HPV"
   - "Human immunodeficiency virus" -> "HIV"
   - "Human herpesvirus"      -> "HHV"
   - "Porcine endogenous retrovirus" -> "PERV"

Usage:
  python3 scripts/shorten_viral_names.py gs://<bucket>/refs/HumanViral_Reference_02-07-2022_modified.renamed_map.tsv
  python3 scripts/shorten_viral_names.py /mnt/disks/staff/refs/HumanViral_Reference_02-07-2022_modified.renamed_map.tsv
"""

import sys
import os
import argparse
import subprocess
import tempfile

REMOVALS = [
    " genomic RNA, complete genome, strain:",
    ", complete genome, strain:",
    " genomic RNA, complete",
    " genomic DNA, partial",
    " genomic DNA.",
    ", complete genome, ",
    ", complete genome,",
    ", complete genome.",
    ", complete genome",
    " complete genome,",
    " complete genome.",
    " complete genome",
    "sense strand",
]

REPLACEMENTS = [
    ("Torque teno mini virus", "TTV mini"),
    ("Torque teno virus", "TTV"),
    ("Adeno-associated virus", "AAV"),
    ("Simian adenovirus", "Simian Ad"),
    ("Human papillomavirus", "HPV"),
    ("Human immunodeficiency virus", "HIV"),
    ("Human herpesvirus", "HHV"),
    ("Porcine endogenous retrovirus", "PERV"),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Shorten viral descriptive names into 3-column .renamed_map_short.tsv file."
    )
    parser.add_argument(
        "input_map",
        help="Path or GCS URI to input rename map TSV (e.g. gs://<bucket>/refs/HumanViral_Reference_02-07-2022_modified.renamed_map.tsv)",
    )
    return parser.parse_args()


def derive_output_path(input_path):
    """Automatically constructs the output path with '_short.tsv' in the exact same directory."""
    if input_path.endswith(".tsv"):
        return input_path[:-4] + "_short.tsv"
    return input_path + "_short.tsv"


def shorten_text(text):
    """Applies rule-based removals and replacements to a string."""
    result = text

    # 1. Apply removals
    for target in REMOVALS:
        result = result.replace(target, "")

    # 2. Apply replacements
    for old_term, new_term in REPLACEMENTS:
        result = result.replace(old_term, new_term)

    # Clean up leftover double spaces or trailing punctuation
    result = " ".join(result.split())
    result = result.strip(",. ")
    return result


def is_gcs_path(path):
    return path.startswith("gs://")


def download_from_gcs(gcs_uri, local_path):
    print(f"[INFO] Downloading input map from GCS: {gcs_uri} ...")
    cmd = f"gsutil -q cp \"{gcs_uri}\" \"{local_path}\""
    subprocess.run(cmd, shell=True, check=True)


def upload_to_gcs(local_path, gcs_uri):
    print(f"[INFO] Uploading shortened map file to GCS: {gcs_uri} ...")
    cmd = f"gsutil -q cp \"{local_path}\" \"{gcs_uri}\""
    subprocess.run(cmd, shell=True, check=True)


def process_map_file(local_input, local_output):
    print(f"[INFO] Processing rename map file into 3-column format...")
    rows_processed = 0

    with open(local_input, "r", encoding="utf-8", errors="replace") as infile, open(
        local_output, "w", encoding="utf-8", newline="\n"
    ) as outfile:
        header_seen = False
        for line in infile:
            line_str = line.rstrip("\r\n")
            if not line_str:
                continue

            parts = line_str.split("\t")
            if not header_seen and ("final_clean_id" in line_str or "original_header" in line_str):
                header_seen = True
                outfile.write("final_clean_id\toriginal_header\tshort_name\n")
                continue

            if len(parts) >= 2:
                clean_id = parts[0]
                orig_header = parts[1]

                # Strip leading accession ID if present at the start of original header
                if orig_header.startswith(clean_id):
                    desc = orig_header[len(clean_id):].strip()
                else:
                    desc = orig_header

                # Apply shortening rules to the description (without accession)
                short_name = shorten_text(desc)
                if not short_name:
                    short_name = clean_id

                outfile.write(f"{clean_id}\t{orig_header}\t{short_name}\n")
                rows_processed += 1
            else:
                outfile.write(line_str + "\n")

    print(f"[SUCCESS] Processed {rows_processed} entries into 3-column format.")


def main():
    args = parse_args()
    input_path = args.input_map
    output_path = derive_output_path(input_path)

    print(f"=== VIROnator Viral Name Shortener ===")
    print(f"Input Map Path:  {input_path}")
    print(f"Output Map Path: {output_path}")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Step 1: Download/resolve input file
        if is_gcs_path(input_path):
            local_in = os.path.join(tmpdir, "input_map.tsv")
            download_from_gcs(input_path, local_in)
        else:
            local_in = input_path

        if not os.path.exists(local_in):
            sys.exit(f"Error: Input file does not exist: {local_in}")

        # Step 2: Setup local output path
        local_out = os.path.join(tmpdir, "output_map.tsv") if is_gcs_path(output_path) else output_path
        os.makedirs(os.path.dirname(os.path.abspath(local_out)), exist_ok=True)

        # Step 3: Shorten names into 3-column TSV
        process_map_file(local_in, local_out)

        # Step 4: Upload to GCS or verify local output
        if is_gcs_path(output_path):
            upload_to_gcs(local_out, output_path)
            print(f"[SUCCESS] Uploaded 3-column shortened map to GCS: {output_path}")
        else:
            print(f"[SUCCESS] Saved 3-column shortened map to local file: {output_path}")


if __name__ == "__main__":
    main()
