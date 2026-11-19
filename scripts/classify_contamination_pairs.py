#!/usr/bin/env python3
"""
scripts/classify_contamination_pairs.py
=======================================
Implements paired-end contamination classification decision tree for VIROnator Step 2.
Classifies each mate (R1 and R2) of a read pair across:
  - VIRUS
  - VECTOR (passes QC: aligned_len >= min_len AND NM <= max_nm)
  - VECTOR_WEAK (fails QC: aligned_len < min_len OR NM > max_nm)
  - MOUSE (passes QC: (aligned_len >= min_len AND NM <= max_nm) OR MAPQ >= min_mapq)
  - MOUSE_WEAK (fails QC)
  - UNMAPPED

Emits:
  - viral_clean.keep: Read IDs to retain for clean viral candidate FASTQs
  - discarded_read_ids.tsv.gz: 2-column gzipped audit table of discarded read IDs and their category
  - distribution_reads_summary_mapq{mapq}_nm{nm}_len{len}.tsv: Summary table of individual read mates
  - distribution_pairs_summary_mapq{mapq}_nm{nm}_len{len}.tsv: Full audit table of read pair combinations
"""

import sys
import re
import os
import gzip
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
    parser.add_argument("--plasmid-contigs", required=True, help="1-column vector/plasmid contig names")
    parser.add_argument("--mouse-contigs", required=True, help="1-column mouse contig names")
    parser.add_argument("--min-len", type=int, default=50, help="Minimum aligned bases (default: 50)")
    parser.add_argument("--max-nm", type=int, default=3, help="Maximum NM edit distance (default: 3, for NM <= 3)")
    parser.add_argument("--min-mapq", type=int, default=20, help="Minimum MAPQ for confident host/mouse (default: 20)")
    parser.add_argument("--keep-plasmid-unmapped", action="store_true", default=False, help="Keep vector-unmapped pairs (default: False)")
    parser.add_argument("--out-clean-keep", required=True, help="Output file for clean viral read IDs")
    parser.add_argument("--out-discarded-ids-gz", required=True, help="Output file for gzipped discarded read IDs (2-column TSV)")
    parser.add_argument("--out-reads-summary", required=True, help="Output file for individual read breakdown TSV")
    parser.add_argument("--out-pairs-summary", required=True, help="Output file for pair combinations breakdown TSV")

    args = parser.parse_args()

    viral_contigs = load_contig_set(args.viral_contigs)
    vector_contigs = load_contig_set(args.plasmid_contigs)
    mouse_contigs = load_contig_set(args.mouse_contigs)

    print(f"[CLASSIFY] Loaded {len(viral_contigs)} Viral, {len(vector_contigs)} Vector, {len(mouse_contigs)} Mouse contigs.")
    print(f"[CLASSIFY] Thresholds: min_len={args.min_len}, max_nm={args.max_nm}, min_mapq={args.min_mapq}")

    # Track alignments per mate: pair_id -> {1: label, 2: label}
    # Labels: 'VIRUS', 'VECTOR', 'VECTOR_WEAK', 'MOUSE', 'MOUSE_WEAK', 'UNMAPPED'
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
            elif rname in vector_contigs:
                aligned_len = parse_cigar_aligned_len(cigar)
                nm = parse_nm_tag(tags)
                if qlen >= args.min_len and aligned_len >= args.min_len and nm != -1 and nm <= args.max_nm:
                    label = "VECTOR"
                else:
                    label = "VECTOR_WEAK"
            elif rname in mouse_contigs:
                aligned_len = parse_cigar_aligned_len(cigar)
                nm = parse_nm_tag(tags)
                if (qlen >= args.min_len and aligned_len >= args.min_len and nm != -1 and nm <= args.max_nm) or mapq >= args.min_mapq:
                    label = "MOUSE"
                else:
                    label = "MOUSE_WEAK"
            else:
                label = "UNMAPPED"

            # Assign label; prioritize high-confidence hits if secondary alignments exist
            prev = mate_labels[pair_id].get(mate)
            priority = {"VIRUS": 5, "VECTOR": 4, "MOUSE": 4, "VECTOR_WEAK": 2, "MOUSE_WEAK": 2, "UNMAPPED": 0}
            if prev is None or priority.get(label, 0) > priority.get(prev, 0):
                mate_labels[pair_id][mate] = label

    # Count individual read mates
    total_pairs = len(mate_labels)
    total_reads = total_pairs * 2
    read_counts = defaultdict(int)

    for pid, mates in mate_labels.items():
        m1 = mates.get(1, "UNMAPPED")
        m2 = mates.get(2, "UNMAPPED")
        read_counts[m1] += 1
        read_counts[m2] += 1

    # Decision tree rules definition
    # Key: (sorted_mate1, sorted_mate2) -> (action, rationale)
    DECISION_RULES = {
        ("VIRUS", "VIRUS"): ("KEEP", "Confirmed viral pair"),
        ("MOUSE", "VIRUS"): ("KEEP", "Host-virus integration junction"),
        ("MOUSE_WEAK", "VIRUS"): ("KEEP", "Potential divergent integration junction"),
        ("UNMAPPED", "VIRUS"): ("KEEP", "Divergent viral strain or novel junction"),
        ("VECTOR_WEAK", "VIRUS"): ("KEEP", "Short micro-homology; viral signal retained for Step 4/6"),
        ("VECTOR", "VIRUS"): ("DISCARD", "Recombinant vector backbone carrying viral promoter/gene"),
        ("VECTOR", "VECTOR"): ("DISCARD", "Pure vector contaminant"),
        ("VECTOR", "VECTOR_WEAK"): ("DISCARD", "Vector contaminant with one degraded mate"),
        ("MOUSE", "VECTOR"): ("DISCARD", "Vector-host cloning artifact"),
        ("MOUSE_WEAK", "VECTOR"): ("DISCARD", "Vector-host cloning artifact"),
        ("UNMAPPED", "VECTOR"): ("KEEP" if args.keep_plasmid_unmapped else "DISCARD", "Vector fragment or cloning boundary"),
        ("MOUSE", "MOUSE"): ("DISCARD", "Host genomic DNA"),
        ("MOUSE", "MOUSE_WEAK"): ("DISCARD", "Host genomic DNA with divergent mate"),
        ("MOUSE", "UNMAPPED"): ("DISCARD", "Host genomic fragment"),
        ("VECTOR_WEAK", "VECTOR_WEAK"): ("DISCARD", "Weak vector noise with no viral evidence"),
        ("MOUSE_WEAK", "VECTOR_WEAK"): ("DISCARD", "Weak vector/host noise with no viral evidence"),
        ("UNMAPPED", "VECTOR_WEAK"): ("DISCARD", "Weak vector noise with no viral evidence"),
        ("MOUSE_WEAK", "MOUSE_WEAK"): ("DISCARD", "Weak host noise with no viral evidence"),
        ("MOUSE_WEAK", "UNMAPPED"): ("DISCARD", "Weak host noise with no viral evidence"),
        ("UNMAPPED", "UNMAPPED"): ("DISCARD", "Uninformative sequencing noise")
    }

    clean_keep = set()
    discarded_records = []
    pair_counts = defaultdict(int)

    for pair_id, mates in mate_labels.items():
        m1 = mates.get(1, "UNMAPPED")
        m2 = mates.get(2, "UNMAPPED")

        sorted_pair = tuple(sorted([m1, m2]))
        action, rationale = DECISION_RULES.get(sorted_pair, ("DISCARD", "Unclassified sequencing noise"))

        pair_counts[sorted_pair] += 1

        if action == "KEEP":
            clean_keep.add(pair_id)
        else:
            cat_str = f"{sorted_pair[0]}-{sorted_pair[1]}"
            discarded_records.append((pair_id, cat_str))

    # 1. Write clean viral IDs
    with open(args.out_clean_keep, 'w') as f_clean:
        for pid in sorted(clean_keep):
            f_clean.write(f"{pid}\n")

    # 2. Write gzipped discarded read IDs audit table
    with gzip.open(args.out_discarded_ids_gz, 'wt') as f_disc:
        f_disc.write("read_id\tpair_category\n")
        for pid, cat in discarded_records:
            f_disc.write(f"{pid}\t{cat}\n")

    # 3. Write individual reads breakdown TSV
    ordered_read_types = ["VIRUS", "VECTOR", "VECTOR_WEAK", "MOUSE", "MOUSE_WEAK", "UNMAPPED"]
    with open(args.out_reads_summary, 'w') as f_reads:
        f_reads.write("mate_type\tcount\tpercentage_of_total\n")
        f_reads.write(f"TOTAL\t{total_reads}\t100.00%\n")
        for rtype in ordered_read_types:
            cnt = read_counts.get(rtype, 0)
            pct = (cnt / total_reads * 100.0) if total_reads > 0 else 0.0
            f_reads.write(f"{rtype}\t{cnt}\t{pct:.2f}%\n")

    # 4. Write read pairs combinations breakdown TSV
    # Define display order
    ordered_pair_types = [
        ("VIRUS", "VIRUS"),
        ("MOUSE", "VIRUS"),
        ("MOUSE_WEAK", "VIRUS"),
        ("UNMAPPED", "VIRUS"),
        ("VECTOR_WEAK", "VIRUS"),
        ("VECTOR", "VIRUS"),
        ("VECTOR", "VECTOR"),
        ("VECTOR", "VECTOR_WEAK"),
        ("MOUSE", "VECTOR"),
        ("MOUSE_WEAK", "VECTOR"),
        ("UNMAPPED", "VECTOR"),
        ("MOUSE", "MOUSE"),
        ("MOUSE", "MOUSE_WEAK"),
        ("MOUSE", "UNMAPPED"),
        ("VECTOR_WEAK", "VECTOR_WEAK"),
        ("MOUSE_WEAK", "VECTOR_WEAK"),
        ("UNMAPPED", "VECTOR_WEAK"),
        ("MOUSE_WEAK", "MOUSE_WEAK"),
        ("MOUSE_WEAK", "UNMAPPED"),
        ("UNMAPPED", "UNMAPPED"),
    ]

    with open(args.out_pairs_summary, 'w') as f_pairs:
        f_pairs.write("R1_mate_type\tR2_mate_type\taction\trationale\tcount\tpercentage_of_total\n")
        f_pairs.write(f"TOTAL\tTOTAL\tTOTAL\tTotal candidate read pairs entering Step 1\t{total_pairs}\t100.00%\n")
        for ptype in ordered_pair_types:
            cnt = pair_counts.get(ptype, 0)
            pct = (cnt / total_pairs * 100.0) if total_pairs > 0 else 0.0
            action, rationale = DECISION_RULES.get(ptype, ("DISCARD", "Other"))
            f_pairs.write(f"{ptype[0]}\t{ptype[1]}\t{action}\t{rationale}\t{cnt}\t{pct:.2f}%\n")

    print(f"\n[CLASSIFY SUMMARY] Total Candidate Pairs Processed: {total_pairs}")
    print(f"[CLASSIFY SUMMARY] Total Clean Viral Pairs Retained : {len(clean_keep)}")
    print(f"[CLASSIFY SUMMARY] Total Discarded Pairs Recorded   : {len(discarded_records)}")
    print(f"[CLASSIFY SUMMARY] Output written:")
    print(f"  - {args.out_clean_keep}")
    print(f"  - {args.out_discarded_ids_gz}")
    print(f"  - {args.out_reads_summary}")
    print(f"  - {args.out_pairs_summary}")

if __name__ == "__main__":
    main()
