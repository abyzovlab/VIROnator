#!/usr/bin/env python3
"""
parse_kraken_report.py

Parses Kraken2 standard hierarchical reports (_kraken2_report.txt) per sample and aggregates
cohort-wide findings into a standardized tabular format directly comparable with VIROnator
and CZ-ID reports.

Features:
1. Reconstructs taxonomic hierarchy (species -> genus -> family -> realm/domain).
2. Filters for viral clades (NCBI TaxID 10239 Viruses) and flags plasmid / synthetic vector hits.
3. Produces per-sample parsed summary TSVs.
4. Aggregates cohort samples into a unified Kraken2 master report.
5. Optionally performs tri-way comparative merging with VIROnator and CZ-ID master reports.
"""

import os
import sys
import argparse
import csv
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Parse and aggregate Kraken2 reports for VIROnator benchmarking.")
    parser.add_argument("--kraken-dir", required=True, help="Directory containing per-sample *_kraken2_report.txt files.")
    parser.add_argument("--out-report", default="cohort_kraken2_master_report.tsv", help="Output path for aggregated Kraken2 report.")
    parser.add_argument("--min-reads", type=int, default=1, help="Minimum cumulative clade reads required to report a taxon (default: 1).")
    parser.add_argument("--target-rank", default="S", help="Target taxonomic rank to report: S (species), G (genus), or all (default: S).")
    parser.add_argument("--filter-viruses-only", action="store_true", default=True, help="Only report taxa within the Viruses clade (taxid 10239) or non-human vectors.")
    parser.add_argument("--merge-vironator-report", default=None, help="Optional path to VIROnator/CZ-ID master report to produce a tri-way concordance table.")
    parser.add_argument("--tri-output", default="cohort_tri_comparison_master.tsv", help="Output path for tri-way comparison TSV if merged.")
    return parser.parse_args()

def parse_single_kraken_report(report_path, min_reads=1, target_rank='S', filter_viruses_only=True):
    """
    Parses a single standard 6-column Kraken2 report.
    Format:
      col 0: percent of fragments (float)
      col 1: cumulative fragments covered by clade (int)
      col 2: fragments assigned directly to taxon (int)
      col 3: rank code (U, R, D, K, P, C, O, F, G, S, etc.)
      col 4: NCBI taxid (str)
      col 5: scientific name (indented with 2 spaces per tree depth level)
    """
    results = []
    lineage_stack = []  # Stack of (indent_level, rank, taxid, name)
    is_viral_branch = False

    with open(report_path, 'r', encoding='utf-8') as f:
        for line in f:
            line_str = line.rstrip('\r\n')
            if not line_str:
                continue
            parts = line_str.split('\t')
            if len(parts) < 6:
                continue

            pct = float(parts[0].strip())
            clade_reads = int(parts[1].strip())
            direct_reads = int(parts[2].strip())
            rank_code = parts[3].strip()
            taxid = parts[4].strip()
            raw_name = parts[5]

            # Calculate indent level (leading spaces indicate tree depth)
            indent_level = len(raw_name) - len(raw_name.lstrip(' '))
            clean_name = raw_name.strip()

            # Maintain lineage stack
            while lineage_stack and lineage_stack[-1][0] >= indent_level:
                lineage_stack.pop()

            current_node = (indent_level, rank_code, taxid, clean_name)
            lineage_stack.append(current_node)

            # Check if this node or any ancestor is Viral (taxid 10239)
            ancestor_taxids = [x[2] for x in lineage_stack]
            ancestor_names = [x[3].lower() for x in lineage_stack]
            
            is_virus = ("10239" in ancestor_taxids or 
                        "viruses" in ancestor_names or 
                        "viroid" in ancestor_names or
                        "plasmid" in clean_name.lower() or
                        "vector" in clean_name.lower())

            # Skip unclassified (taxid 0) and host (human taxid 9606)
            if taxid == "0" or "9606" in ancestor_taxids:
                continue

            if filter_viruses_only and not is_virus:
                continue

            if clade_reads < min_reads:
                continue

            # Check target rank match (e.g. 'S' or 'S1' for species / strain)
            rank_match = False
            if target_rank == 'all':
                rank_match = True
            elif target_rank == 'S' and (rank_code == 'S' or rank_code.startswith('S')):
                rank_match = True
            elif rank_code == target_rank:
                rank_match = True

            if rank_match:
                # Extract genus and family from lineage stack
                genus_taxid, genus_name = "NA", "NA"
                family_taxid, family_name = "NA", "NA"
                for node in lineage_stack:
                    if node[1] == 'G':
                        genus_taxid, genus_name = node[2], node[3]
                    elif node[1] == 'F':
                        family_taxid, family_name = node[2], node[3]

                results.append({
                    'species_taxid': taxid,
                    'species_name': clean_name,
                    'rank': rank_code,
                    'kraken_clade_reads': clade_reads,
                    'kraken_direct_reads': direct_reads,
                    'kraken_percent': pct,
                    'genus_taxid': genus_taxid,
                    'genus_name': genus_name,
                    'family_taxid': family_taxid,
                    'family_name': family_name
                })

    return results

def main():
    args = parse_args()
    kraken_dir = args.kraken_dir
    if not os.path.exists(kraken_dir):
        print(f"[ERROR] Directory not found: {kraken_dir}", file=sys.stderr)
        sys.exit(1)

    # Find all *_kraken2_report.txt files recursively
    report_files = []
    for root, _, files in os.walk(kraken_dir):
        for f in files:
            if f.endswith("_kraken2_report.txt") or f.endswith(".report"):
                report_files.append(os.path.join(root, f))

    if not report_files:
        print(f"[WARNING] No *_kraken2_report.txt files found in {kraken_dir}!", file=sys.stderr)
        sys.exit(0)

    print(f"[INFO] Discovered {len(report_files)} Kraken2 report file(s). Parsing...")

    cohort_rows = []
    sample_to_kraken = defaultdict(list)

    for r_path in sorted(report_files):
        # Extract sample_id from filename (e.g. s_3595-ST-0002_kraken2_report.txt -> s_3595-ST-0002)
        base = os.path.basename(r_path)
        sample_id = base.replace("_kraken2_report.txt", "").replace(".report", "")

        parsed_taxa = parse_single_kraken_report(
            r_path,
            min_reads=args.min_reads,
            target_rank=args.target_rank,
            filter_viruses_only=args.filter_viruses_only
        )

        for taxon in parsed_taxa:
            row = {'sample_id': sample_id}
            row.update(taxon)
            row['source_report'] = r_path
            cohort_rows.append(row)
            sample_to_kraken[sample_id].append(row)

    # Write aggregated master Kraken2 report
    out_headers = [
        'sample_id', 'species_taxid', 'species_name', 'rank',
        'kraken_clade_reads', 'kraken_direct_reads', 'kraken_percent',
        'genus_taxid', 'genus_name', 'family_taxid', 'family_name', 'source_report'
    ]

    with open(args.out_report, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=out_headers, delimiter='\t')
        writer.writeheader()
        for r in cohort_rows:
            writer.writerow(r)

    print(f"[SUCCESS] Aggregated {len(cohort_rows)} viral detections across {len(sample_to_kraken)} sample(s).")
    print(f"[OUTPUT] Master report saved to: {args.out_report}")

    # Optional tri-way comparison merge with VIROnator / CZ-ID master report
    if args.merge_vironator_report and os.path.exists(args.merge_vironator_report):
        print(f"[INFO] Merging with VIROnator/CZ-ID master report: {args.merge_vironator_report}...")
        merge_tri_comparison(args.merge_vironator_report, sample_to_kraken, args.tri_output)

def merge_tri_comparison(vironator_master_path, sample_to_kraken, output_path):
    """
    Merges VIROnator, CZ-ID, and Kraken2 calls side-by-side by matching sample_id and species_taxid.
    """
    with open(vironator_master_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        v_rows = list(reader)

    tri_rows = []
    all_samples = set([r['sample_id'] for r in v_rows] + list(sample_to_kraken.keys()))

    # Index VIROnator rows by (sample_id, species_taxid)
    v_by_sample_sp = defaultdict(list)
    for r in v_rows:
        sid = r['sample_id']
        sp_id = r.get('species_taxid', '').strip()
        v_by_sample_sp[(sid, sp_id)].append(r)

    processed_pairs = set()

    for r in v_rows:
        sid = r['sample_id']
        sp_id = r.get('species_taxid', '').strip()
        pair_key = (sid, sp_id)
        processed_pairs.add(pair_key)

        # Look for matching taxon in Kraken2 results for this sample
        matching_kraken = [k for k in sample_to_kraken.get(sid, []) if k['species_taxid'] == sp_id]
        if matching_kraken:
            k_top = sorted(matching_kraken, key=lambda x: x['kraken_clade_reads'], reverse=True)[0]
            k_clade = k_top['kraken_clade_reads']
            k_direct = k_top['kraken_direct_reads']
            k_pct = k_top['kraken_percent']
        else:
            k_clade, k_direct, k_pct = 0, 0, 0.0

        tri_rows.append({
            'sample_id': sid,
            'species_taxid': sp_id,
            'species_name': r.get('species_name', ''),
            'clinician_reported': r.get('clinician_reported', 'NA'),
            'vironator_accession': r.get('virus_accession', 'NA'),
            'vironator_reads': r.get('virus_mapped_reads', 0),
            'vironator_copy_number': r.get('viral_copy_number', 0.0),
            'vironator_physical_cov': r.get('physical_coverage', 0.0),
            'czid_virus_name': r.get('call_virus_name', 'NA'),
            'czid_species_taxid': r.get('call_species_taxid', 'NA'),
            'czid_reads': r.get('call_species_total_reads', 0),
            'kraken_clade_reads': k_clade,
            'kraken_direct_reads': k_direct,
            'kraken_percent': k_pct,
            'concordance_all_three': 'YES' if (float(r.get('virus_mapped_reads', 0)) > 0 and float(r.get('call_species_total_reads', 0)) > 0 and k_clade > 0) else 'NO'
        })

    # Add Kraken2 calls not observed in VIROnator
    for sid, k_list in sample_to_kraken.items():
        for k in k_list:
            sp_id = k['species_taxid']
            if (sid, sp_id) not in processed_pairs:
                tri_rows.append({
                    'sample_id': sid,
                    'species_taxid': sp_id,
                    'species_name': k['species_name'],
                    'clinician_reported': 'NA',
                    'vironator_accession': 'None',
                    'vironator_reads': 0,
                    'vironator_copy_number': 0.0,
                    'vironator_physical_cov': 0.0,
                    'czid_virus_name': 'None',
                    'czid_species_taxid': 'NA',
                    'czid_reads': 0,
                    'kraken_clade_reads': k['kraken_clade_reads'],
                    'kraken_direct_reads': k['kraken_direct_reads'],
                    'kraken_percent': k['kraken_percent'],
                    'concordance_all_three': 'NO'
                })

    headers = [
        'sample_id', 'species_taxid', 'species_name', 'clinician_reported',
        'vironator_accession', 'vironator_reads', 'vironator_copy_number', 'vironator_physical_cov',
        'czid_virus_name', 'czid_species_taxid', 'czid_reads',
        'kraken_clade_reads', 'kraken_direct_reads', 'kraken_percent',
        'concordance_all_three'
    ]

    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers, delimiter='\t')
        writer.writeheader()
        for tr in tri_rows:
            writer.writerow(tr)

    print(f"[SUCCESS] Tri-way comparative report saved to: {output_path}")

if __name__ == '__main__':
    main()
