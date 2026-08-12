#!/usr/bin/env python3
"""
generate_report.py
Consolidated viral reporting script for VIROnator pipeline.

Generates a 14-column TSV report per sample evaluating up to 4 CRAM output files:
- exogeneSR_viral_clean.sorted.cram
- exogeneSR_viral_clean.sorted.flags.cram
- exogeneSR_viral_raw.sorted.cram
- exogeneSR_viral_raw.sorted.flags.cram

Columns:
1. Phase
2. Project (or 'base' if empty)
3. Source_File
4. Sample_ID
5. Virus_Accession
6. Virus_Length
7. Virus_Mapped_Reads
8. Normalized_Coverage
9. Physical_Coverage
10. Human_Genome_Size
11. Sample_Read_Depth
12. Viral_Copy_Number
13. Specimen
14. Virus_Name_Sanitized
"""

import argparse
import os
import re
import subprocess
import sys


def parse_args():
    parser = argparse.ArgumentParser(description="Generate consolidated viral report.")
    parser.add_argument("--sample-id", required=True, help="Sample ID (e.g. SSC01120)")
    parser.add_argument("--phase", required=True, help="Pipeline phase (e.g. 2)")
    parser.add_argument("--project", default="", help="Project name (defaults to 'base' if empty)")
    parser.add_argument("--vironator-dir", required=True, help="Path to sample vironator output folder")
    parser.add_argument("--out-file", required=True, help="Path to output TSV report file")
    parser.add_argument("--combined-ref", required=True, help="Path to combined reference FASTA")
    parser.add_argument("--human-ref-fai", required=True, help="Path to human reference FASTA .fai index")
    parser.add_argument("--rename-map", required=True, help="Path to viral rename_map.tsv file")
    parser.add_argument("--viral-bed", required=True, help="Path to viral BED file")
    parser.add_argument("--metadata", required=True, help="Path to SSC_sample_metadata.tsv")
    return parser.parse_args()


def load_viral_lengths(bed_filepath, fai_filepath=""):
    """Loads viral contig length mapping from BED file (column 3) or FASTA .fai index."""
    lengths = {}
    if os.path.exists(bed_filepath):
        with open(bed_filepath, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    contig_id = parts[0]
                    try:
                        lengths[contig_id] = int(parts[2])
                    except ValueError:
                        continue

    if not lengths and fai_filepath and os.path.exists(fai_filepath):
        with open(fai_filepath, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    contig_id = parts[0]
                    try:
                        lengths[contig_id] = int(parts[1])
                    except ValueError:
                        continue
    return lengths


def load_viral_names(rename_map_filepath):
    """Loads viral contig descriptive name mapping from rename_map.tsv (final_clean_id -> original_header)."""
    names = {}
    if not os.path.exists(rename_map_filepath):
        return names
    with open(rename_map_filepath, "r") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header is None and ("final_clean_id" in line or "original_header" in line):
                header = [p.lower() for p in parts]
                continue

            if len(parts) >= 2:
                clean_id = parts[0]
                orig_header = parts[1].strip()
                
                # Strip accession prefix if present at start of header
                if orig_header.startswith(clean_id):
                    desc = orig_header[len(clean_id):].strip()
                else:
                    desc = orig_header
                
                # Sanitize description: replace non-alphanumeric characters with underscores
                import re
                sanitized = re.sub(r"[^A-Za-z0-9]+", "_", desc).strip("_")
                names[clean_id] = sanitized if sanitized else clean_id
    return names


def calculate_human_genome_size(fai_filepath):
    """Calculates total human genome size by summing contig lengths in .fai file."""
    if not os.path.exists(fai_filepath):
        return 3099734149  # Default GRCh38 human genome size in bp
    total_size = 0
    with open(fai_filepath, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    total_size += int(parts[1])
                except ValueError:
                    continue
    return total_size if total_size > 0 else 3099734149


def load_metadata(filepath, sample_id, phase, project):
    """
    Loads sample read depth and specimen from sample_metadata.tsv.
    Expected columns: sample, coverage, specimen, phase, project
    """
    default_depth = 30.0
    default_specimen = "Unknown"
    
    project_key = project if project.strip() else "base"
    clean_sample = sample_id.replace(".sorted", "").replace("Sample_", "")
    
    # Format phase string (e.g. phase1 or 1 as passed)
    target_phase = str(phase).strip()
    if not target_phase.startswith("phase") and target_phase.isdigit():
        target_phase_alt = f"phase{target_phase}"
    else:
        target_phase_alt = target_phase

    if not os.path.exists(filepath):
        return default_depth, default_specimen

    with open(filepath, "r") as f:
        header = None
        for line in f:
            parts = line.strip().split("\t")
            if not parts or not line.strip():
                continue
            if header is None:
                header = [p.lower() for p in parts]
                continue
            
            row = dict(zip(header, parts))
            row_sample = row.get("sample", "").replace(".sorted", "").replace("Sample_", "")
            row_phase = str(row.get("phase", "")).strip()
            row_project = str(row.get("project", "")).strip()
            
            if row_sample == clean_sample and (row_phase == target_phase or row_phase == target_phase_alt) and row_project == project_key:
                try:
                    depth = float(row.get("coverage", 30.0))
                except ValueError:
                    depth = 30.0
                specimen = row.get("specimen", "Unknown")
                return depth, specimen

    return default_depth, default_specimen


def get_mapped_reads_per_virus(cram_path, combined_ref):
    """
    Returns a dict of {viral_contig: read_pair_count} from CRAM file,
    strictly counting read pairs where BOTH mates mapped to the exact same viral contig.
    """
    cmd = f"samtools view -f 2 -F 3844 -T {combined_ref} {cram_path} 2>/dev/null | awk '($7 == \"=\" || $7 == $3) && $3 != \"*\" {{print $3}}' | sort | uniq -c"
    try:
        res = subprocess.check_output(cmd, shell=True, text=True)
        counts = {}
        for line in res.strip().split("\n"):
            parts = line.strip().split()
            if len(parts) == 2:
                raw_cnt = int(parts[0])
                contig = parts[1]
                # Each read pair has 2 read records in SAM (R1 and R2), so pair count = raw_cnt / 2
                pair_cnt = max(1, raw_cnt // 2)
                if contig != "*":
                    counts[contig] = pair_cnt
        return counts
    except Exception:
        return {}


def calculate_physical_coverage(cram_path, combined_ref, virus, virus_size):
    """Calculates physical coverage as a percentage (0.00 to 100.00) without % sign."""
    if virus_size <= 0:
        return "00.00"
    # Fast region-bounded mpileup directly targeting target viral contig -r {virus}
    cmd = f"samtools mpileup -r '{virus}' -f {combined_ref} {cram_path} 2>/dev/null | cut -f2,4"
    try:
        res = subprocess.check_output(cmd, shell=True, text=True)
    except subprocess.CalledProcessError:
        return "00.00"

    covered_bases = 0
    for line in res.strip().split("\n"):
        parts = line.strip().split()
        if len(parts) >= 2:
            try:
                depth = int(parts[1])
                if depth > 0:
                    covered_bases += 1
            except ValueError:
                continue

    pct = (covered_bases / float(virus_size)) * 100.0
    pct = min(pct, 100.0)
    return f"{pct:05.2f}"


def main():
    args = parse_args()
    
    project_label = args.project.strip() if args.project.strip() else "base"
    phase_label = f"phase{args.phase}" if not str(args.phase).startswith("phase") else str(args.phase)
    viral_lengths = load_viral_lengths(args.viral_bed)
    viral_names = load_viral_names(args.rename_map)
    human_genome_size = calculate_human_genome_size(args.human_ref_fai)
    read_depth, specimen = load_metadata(args.metadata, args.sample_id, args.phase, args.project)

    target_files = [
        ("exogeneSR_viral_clean.sorted.cram", os.path.join(args.vironator_dir, "exogeneSR_viral_clean.sorted.cram")),
        ("exogeneSR_viral_clean.sorted.flags.cram", os.path.join(args.vironator_dir, "exogeneSR_viral_clean.sorted.flags.cram")),
        ("exogeneSR_viral_raw.sorted.cram", os.path.join(args.vironator_dir, "exogeneSR_viral_raw.sorted.cram")),
        ("exogeneSR_viral_raw.sorted.flags.cram", os.path.join(args.vironator_dir, "exogeneSR_viral_raw.sorted.flags.cram")),
    ]

    os.makedirs(os.path.dirname(args.out_file), exist_ok=True)

    header = [
        "Sample_ID",
        "Virus_Accession",
        "Virus_Length",
        "Virus_Mapped_Reads",
        "Normalized_Coverage",
        "Physical_Coverage",
        "Human_Genome_Size",
        "Sample_Read_Depth",
        "Viral_Copy_Number",
        "Virus_Name_Sanitized",
        "Specimen",
        "Phase",
        "Project",
        "Source_File",
    ]

    rows = []
    for fname, fpath in target_files:
        if not os.path.exists(fpath):
            continue

        mapped_counts = get_mapped_reads_per_virus(fpath, args.combined_ref)
        positive_hits = {acc: cnt for acc, cnt in mapped_counts.items() if cnt > 0}

        if not positive_hits:
            # Baseline negative record when no viral reads are detected in CRAM file
            row = [
                str(args.sample_id),
                "None",
                "0",
                "0",
                "0.000000",
                "00.00",
                str(human_genome_size),
                f"{read_depth:.2f}",
                "0.000000",
                "None",
                str(specimen),
                phase_label,
                project_label,
                fname,
            ]
            rows.append(row)
        else:
            for virus_accession, read_count in positive_hits.items():
                virus_length = viral_lengths.get(virus_accession, 10000)
                norm_cov = (read_count * 150.0) / float(virus_length)
                phys_cov = calculate_physical_coverage(fpath, args.combined_ref, virus_accession, virus_length)
                
                denom = (read_depth / 2.0) if read_depth > 0 else 15.0
                copy_number = norm_cov * (1.0 / denom)
                virus_name = viral_names.get(virus_accession, virus_accession)

                row = [
                    str(args.sample_id),
                    str(virus_accession),
                    str(virus_length),
                    str(read_count),
                    f"{norm_cov:.6f}",
                    str(phys_cov),
                    str(human_genome_size),
                    f"{read_depth:.2f}",
                    f"{copy_number:.6f}",
                    str(virus_name),
                    str(specimen),
                    phase_label,
                    project_label,
                    fname,
                ]
                rows.append(row)

    with open(args.out_file, "w") as out:
        out.write("\t".join(header) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print(f"Report generated successfully: {args.out_file} ({len(rows)} entries)")


if __name__ == "__main__":
    main()
