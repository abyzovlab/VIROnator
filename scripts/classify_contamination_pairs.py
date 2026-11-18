#!/usr/bin/env python3
"""
scripts/classify_contamination_pairs.py
=======================================
Implements paired-end classification decision tree for VIROnator Step 2.
Classifies each mate (R1 and R2) of a read pair across:
  - Viral contigs
  - Plasmid contigs (aligned_len >= min_len, NM <= max_nm)
  - Mouse contigs (aligned_len >= min_len, NM <= max_nm, or MAPQ >= 20)
  - Unmapped / Sub-threshold

Emits:
  - viral_clean.keep: Read IDs to retain for clean viral candidate FASTQs
  - plasmid_reads_se.ids: Read IDs blacklisted as plasmid
  - mouse_reads_se.ids: Read IDs blacklisted as mouse
  - pair_classification_summary.tsv: Diagnostic audit of pair counts per category
"""

import sys
import re
import argparse
from collections import defaultdict

def parse_cigar_aligned_len(cigar):
    """Calculates total aligned/matched length (M, =, X) from CIGAR string."""
    if not cigar or cigar == "*":
        return 0
    total = 0
    for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in ("M", "=", "X"):
            total += int(length)
    return total

def parse_nm_tag(tags):
    """Extracts NM edit distance integer from SAM optional tags."""
    for tag in tags:
        if tag.startswith("NM:i:"):
            try:
                return int(tag[5:])
            except ValueError:
                return -1
    return -1

def load_contig_set(filepath):
    """Loads 1-column list of contig identifiers."""
    cset = set()
    if filepath and filepath != "none":
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    c = line.strip().split()[0] if line.strip() else ""
                    if c:
                        cset.add(c)
        except FileNotFoundError:
            pass
    return cset

def main():
    parser = argparse.ArgumentParser(description="Paired-end contamination classification decision tree.")
    parser.add_argument("--reads-temp", required=True, help="SAM output from Step 1 single-end screening")
    parser.add_argument("--viral-contigs", required=True, help="1-column viral contig names")
    parser.add_argument("--plasmid-contigs", required=True, help="1-column plasmid contig names")
    parser.add_argument("--mouse-contigs", required=True, help="1-column mouse contig names")
    parser.add_argument("--min-len", type=int, default=50, help="Minimum aligned bases (default: 50)")
    parser.add_argument("--max-nm", type=int, default=3, help="Maximum NM edit distance (default: 3, for NM < 4)")
    parser.add_argument("--keep-plasmid-unmapped", action="store_true", default=True, help="Keep plasmid-unmapped pairs (default: True)")
    parser.add_argument("--out-clean-keep", required=True, help="Output file for clean viral read IDs")
    parser.add_argument("--out-plasmid-ids", required=True, help="Output file for plasmid read IDs")
    parser.add_argument("--out-mouse-ids", required=True, help="Output file for mouse read IDs")
    parser.add_argument("--out-summary", required=False, help="Diagnostic summary file")

    args = parser.parse_args()

    viral_contigs = load_contig_set(args.viral_contigs)
    plasmid_contigs = load_contig_set(args.plasmid_contigs)
    mouse_contigs = load_contig_set(args.mouse_contigs)

    print(f"[CLASSIFY] Loaded {len(viral_contigs)} Viral, {len(plasmid_contigs)} Plasmid, {len(mouse_contigs)} Mouse contigs.")

    # Track alignments per mate: pair_id -> {1: label, 2: label}
    # Labels: 'VIRUS', 'PLASMID', 'MOUSE', 'WEAK_PLASMID', 'WEAK_MOUSE', 'UNMAPPED'
    mate_labels = defaultdict(dict)

    with open(args.reads_temp, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            qname = parts[0]
            flag = int(parts[1])
            rname = parts[2]
            mapq = int(parts[4])
            cigar = parts[5]
            qlen = len(parts[9])
            tags = parts[11:]

            # Determine pair ID and mate number (1 or 2)
            if qname.endswith('/1'):
                pair_id = qname[:-2]
                mate = 1
            elif qname.endswith('/2'):
                pair_id = qname[:-2]
                mate = 2
            elif flag & 64:
                pair_id = qname
                mate = 1
            elif flag & 128:
                pair_id = qname
                mate = 2
            else:
                pair_id = qname
                mate = 1

            # Check mapping status (bit 4: unmapped)
            is_unmapped = (flag & 4) != 0 or rname == "*"

            if is_unmapped:
                label = "UNMAPPED"
            elif rname in viral_contigs:
                label = "VIRUS"
            elif rname in plasmid_contigs:
                aligned_len = parse_cigar_aligned_len(cigar)
                nm = parse_nm_tag(tags)
                if qlen >= args.min_len and aligned_len >= args.min_len and nm != -1 and nm <= args.max_nm:
                    label = "PLASMID"
                else:
                    label = "WEAK_PLASMID"
            elif rname in mouse_contigs:
                aligned_len = parse_cigar_aligned_len(cigar)
                nm = parse_nm_tag(tags)
                if (qlen >= args.min_len and aligned_len >= args.min_len and nm != -1 and nm <= args.max_nm) or mapq >= 20:
                    label = "MOUSE"
                else:
                    label = "WEAK_MOUSE"
            else:
                # Other contig (e.g. human / decoy)
                label = "OTHER"

            # Assign label; prioritize high-confidence hits if secondary alignments exist
            prev = mate_labels[pair_id].get(mate)
            priority = {"VIRUS": 5, "PLASMID": 4, "MOUSE": 4, "WEAK_PLASMID": 2, "WEAK_MOUSE": 2, "OTHER": 1, "UNMAPPED": 0}
            if prev is None or priority.get(label, 0) > priority.get(prev, 0):
                mate_labels[pair_id][mate] = label

    # Decision Tree Evaluation per pair
    clean_keep = set()
    plasmid_ids = set()
    mouse_ids = set()
    category_counts = defaultdict(int)

    for pair_id, mates in mate_labels.items():
        m1 = mates.get(1, "UNMAPPED")
        m2 = mates.get(2, "UNMAPPED")

        # Sort labels for symmetrical decision rules (e.g. VIRUS-MOUSE == MOUSE-VIRUS)
        pair_type = tuple(sorted([m1, m2]))

        # --- Decision Tree Rules ---
        if "VIRUS" in pair_type:
            if pair_type == ("VIRUS", "VIRUS"):
                action = "KEEP"
                category = "virus-virus"
            elif pair_type == ("MOUSE", "VIRUS"):
                action = "KEEP"
                category = "virus-mouse"
            elif pair_type == ("UNMAPPED", "VIRUS"):
                action = "KEEP"
                category = "virus-unmapped"
            elif pair_type == ("PLASMID", "VIRUS"):
                action = "REMOVE"
                category = "virus-plasmid"
                plasmid_ids.add(pair_id)
            elif pair_type == ("VIRUS", "WEAK_PLASMID"):
                action = "KEEP"
                category = "virus-plasmid_weak"
            elif pair_type == ("VIRUS", "WEAK_MOUSE"):
                action = "KEEP"
                category = "virus-mouse_weak"
            else:
                action = "KEEP"
                category = "virus-other"

        elif "PLASMID" in pair_type:
            if pair_type == ("PLASMID", "PLASMID"):
                action = "REMOVE"
                category = "plasmid-plasmid"
                plasmid_ids.add(pair_id)
            elif "MOUSE" in pair_type or "WEAK_MOUSE" in pair_type:
                action = "REMOVE"
                category = "plasmid-mouse"
                plasmid_ids.add(pair_id)
            elif pair_type == ("PLASMID", "WEAK_PLASMID"):
                # Pure plasmid contaminant where mate had slightly more mismatches/clipping
                action = "REMOVE"
                category = "plasmid-plasmid_weak"
                plasmid_ids.add(pair_id)
            elif "UNMAPPED" in pair_type:
                if args.keep_plasmid_unmapped:
                    action = "KEEP"
                    category = "plasmid-unmapped"
                else:
                    action = "REMOVE"
                    category = "plasmid-unmapped"
                    plasmid_ids.add(pair_id)
            else:
                action = "REMOVE"
                category = "plasmid-other"
                plasmid_ids.add(pair_id)

        elif "MOUSE" in pair_type:
            action = "REMOVE"
            category = "mouse-mouse" if pair_type == ("MOUSE", "MOUSE") else "mouse-unmapped"
            mouse_ids.add(pair_id)

        elif "WEAK_PLASMID" in pair_type:
            action = "REMOVE"
            category = "plasmid_weak-other"
            plasmid_ids.add(pair_id)

        elif "WEAK_MOUSE" in pair_type:
            action = "REMOVE"
            category = "mouse_weak-other"
            mouse_ids.add(pair_id)

        else:
            action = "REMOVE"
            category = "unmapped-unmapped"

        category_counts[category] += 1
        if action == "KEEP":
            clean_keep.add(pair_id)

    # Write output files
    with open(args.out_clean_keep, 'w') as f_clean:
        for pid in sorted(clean_keep):
            f_clean.write(f"{pid}\n")

    with open(args.out_plasmid_ids, 'w') as f_plas:
        for pid in sorted(plasmid_ids):
            f_plas.write(f"{pid}\n")

    with open(args.out_mouse_ids, 'w') as f_mouse:
        for pid in sorted(mouse_ids):
            f_mouse.write(f"{pid}\n")

    print(f"\n[CLASSIFY SUMMARY] Total Candidate Pairs Processed: {len(mate_labels)}")
    print(f"[CLASSIFY SUMMARY] Total Clean Viral Pairs Retained : {len(clean_keep)}")
    print(f"[CLASSIFY SUMMARY] Total Plasmid Blacklisted Pairs    : {len(plasmid_ids)}")
    print(f"[CLASSIFY SUMMARY] Total Mouse Blacklisted Pairs      : {len(mouse_ids)}")
    print("\n--- Detailed Category Breakdown ---")
    for cat, count in sorted(category_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {cat:<26}: {count}")

    if args.out_summary:
        with open(args.out_summary, 'w') as f_sum:
            f_sum.write("category\tcount\n")
            for cat, count in sorted(category_counts.items()):
                f_sum.write(f"{cat}\t{count}\n")

if __name__ == "__main__":
    main()
