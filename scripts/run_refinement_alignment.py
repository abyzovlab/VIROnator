#!/usr/bin/env python3
"""
run_refinement_alignment.py
Per-Sample Refinement Alignment, Classification & Reporting Script for VIROnator

Implements the Hierarchical Two-Tier Classification Cascade:
1. Tier 1 (Strict Unique Check):
   - MAPQ >= 20, F 12 == 0 (both mapped), F 2048 == 0 (non-supp), CIGAR 150M/151M, 150 <= |TLEN| <= 1500
   - Uniquely maps to ONE RefSeq accession.
   - Classified as REFERENCE_UNIQUE. (STOP! Do NOT check alternative/secondary alignments or species-supportive!).

2. Tier 2 (Species-Supportive Check):
   - Evaluates ONLY remaining reads (failed Tier 1 strict criteria or multi-mapped).
   - Relaxed MAPQ (MAPQ >= 5), F 12 == 0, F 2048 == 0, 150 <= |TLEN| <= 1500.
   - Tests if all target accessions belong to the SAME species_taxid (zero host/decoy hits).
   - Classified as SPECIES_SUPPORTIVE.

Outputs:
- CRAM: exogeneSR_viral_refinement.cram
- TSV: <phase>_<project>_<sample>_refined_report.tsv (11 Columns)
"""

import os
import sys
import collections
import argparse
import subprocess
import pysam


def parse_args():
    parser = argparse.ArgumentParser(description="Per-sample Bowtie2 alignment, classification and report generation.")
    parser.add_argument("--sample-id", required=True, help="Sample ID")
    parser.add_argument("--phase", required=True, help="Phase number or tag")
    parser.add_argument("--project", default="base", help="Project name (default: base)")
    parser.add_argument("--vironator-dir", required=True, help="Path to sample VIROnator output directory")
    parser.add_argument("--out-dir", required=True, help="Refinement output directory for sample")
    parser.add_argument("--combined-ref", required=True, help="Path to date-stamped combined reference FASTA")
    parser.add_argument("--taxonomy-index", required=True, help="Path to viral_reference_taxonomy_index.tsv")
    parser.add_argument("--bowtie2-bin", default="bowtie2", help="Path to bowtie2 binary")
    parser.add_argument("--max-multimaps", type=int, default=10, help="-k flag for Bowtie2")
    parser.add_argument("--strategy", default="clean", choices=["clean", "raw"], help="Candidate read strategy source")
    parser.add_argument("--human-coverage", type=float, default=30.0, help="Sample human genome read depth coverage")
    return parser.parse_args()


def load_taxonomy_mapping(tax_index_path):
    """Loads accession -> {species_taxid, species_name} from taxonomy index TSV."""
    tax_map = {}
    if not os.path.exists(tax_index_path):
        print(f"[WARNING] Taxonomy index not found: {tax_index_path}")
        return tax_map

    with open(tax_index_path, "r", encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header is None and "accession" in line:
                header = [p.lower() for p in parts]
                continue

            if len(parts) >= 6:
                acc = parts[0]
                acc_ver = parts[1]
                
                sp_tax_idx = header.index("species_taxid") if header and "species_taxid" in header else 4
                sp_nm_idx = header.index("species_name") if header and "species_name" in header else 5

                sp_id = parts[sp_tax_idx] if len(parts) > sp_tax_idx else "NO_SPECIES_RANK"
                sp_nm = parts[sp_nm_idx] if len(parts) > sp_nm_idx else "NO_SPECIES_NAME"

                tax_map[acc] = (sp_id, sp_nm)
                tax_map[acc_ver] = (sp_id, sp_nm)
    return tax_map


def load_reference_lengths(ref_fasta_path):
    """Loads contig lengths from reference FASTA or .fai index."""
    lengths = {}
    fai_path = ref_fasta_path + ".fai"
    if os.path.exists(fai_path):
        with open(fai_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    lengths[parts[0]] = int(parts[1])
    return lengths


def find_candidate_fastqs(vironator_dir, strategy="clean"):
    """Locates candidate FASTQ files (e.g. viral_clean_1.fq, viral_clean_2.fq)."""
    fq1 = os.path.join(vironator_dir, f"viral_{strategy}_1.fq")
    fq2 = os.path.join(vironator_dir, f"viral_{strategy}_2.fq")
    if os.path.exists(fq1) and os.path.exists(fq2) and os.path.getsize(fq1) > 0:
        return fq1, fq2
    return None, None


def main():
    args = parse_args()

    phase_str = f"phase{args.phase}" if not str(args.phase).startswith("phase") else str(args.phase)
    proj_str = str(args.project) if args.project else "base"
    sample_id = str(args.sample_id)

    os.makedirs(args.out_dir, exist_ok=True)
    out_cram = os.path.join(args.out_dir, "exogeneSR_viral_refinement.cram")
    out_report_tsv = os.path.join(args.out_dir, f"{phase_str}_{proj_str}_{sample_id}_refined_report.tsv")

    tax_map = load_taxonomy_mapping(args.taxonomy_index)
    ref_lengths = load_reference_lengths(args.combined_ref)

    fq1, fq2 = find_candidate_fastqs(args.vironator_dir, args.strategy)
    
    if not fq1 or not fq2:
        print(f"[WARNING] Candidate FASTQ files for strategy '{args.strategy}' not found in {args.vironator_dir}. Writing empty report.")
        with open(out_report_tsv, "w", encoding="utf-8") as rf:
            rf.write("Sample_ID\tPhase\tProject\tVirus_Accession\tSpecies_TaxID\tSpecies_Name\tReference_Unique_Reads\tSpecies_Supportive_Reads\tTotal_Refined_Reads\tRefined_Copy_Number\tClassification_Status\n")
        return

    # Run Bowtie2 Multimapping Alignment
    print(f"[ALIGN] Running Bowtie2 multimapping alignment for {sample_id} ...", flush=True)
    bt2_prefix = args.combined_ref
    cmd_align = (
        f"\"{args.bowtie2_bin}\" -x \"{bt2_prefix}\" -1 \"{fq1}\" -2 \"{fq2}\" "
        f"-k {args.max_multimaps} --very-sensitive --no-mixed --no-discordant -p 8 2> \"{args.out_dir}/bowtie2_refinement.log\" "
        f"| samtools view -C -T \"{args.combined_ref}\" - "
        f"| samtools sort -@ 4 -o \"{out_cram}\" -"
    )
    subprocess.run(cmd_align, shell=True, check=True)
    subprocess.run(f"samtools index \"{out_cram}\"", shell=True, check=True)

    print(f"[CLASSIFY] Running Hierarchical Two-Tier Classification Cascade on {out_cram} ...", flush=True)

    # Open CRAM alignment for read classification parsing
    bamfile = pysam.AlignmentFile(out_cram, "rc", reference_filename=args.combined_ref)
    
    # Group read pairs by qname
    read_pairs = collections.defaultdict(list)
    for read in bamfile.fetch(until_eof=True):
        if not read.is_unmapped:
            read_pairs[read.query_name].append(read)

    unique_counts = collections.Counter()      # acc -> unique read count
    supportive_counts = collections.Counter()  # species_taxid -> supportive read count
    species_acc_map = {}                       # species_taxid -> set of accessions

    decoy_host_count = 0
    ambiguous_count = 0

    for qname, alignments in read_pairs.items():
        # Separate R1 and R2 alignments
        r1_alns = [r for r in alignments if r.is_read1]
        r2_alns = [r for r in alignments if r.is_read2]

        if not r1_alns or not r2_alns:
            continue

        # Helper checks for Tier 1 Strict Primary Criteria
        def is_tier1_strict(r):
            if r.mapping_quality < 20 or r.is_secondary or r.is_supplementary:
                return False
            cigar_str = r.cigarstring
            if cigar_str not in ("150M", "151M"):
                return False
            tlen = abs(r.template_length)
            if tlen < 150 or tlen > 1500:
                return False
            return True

        # Check Tier 1 Primary Strict Unique
        tier1_r1 = [r for r in r1_alns if is_tier1_strict(r)]
        tier1_r2 = [r for r in r2_alns if is_tier1_strict(r)]

        is_tier1_pass = False
        if len(tier1_r1) == 1 and len(tier1_r2) == 1:
            ref1 = bamfile.get_reference_name(tier1_r1[0].reference_id)
            ref2 = bamfile.get_reference_name(tier1_r2[0].reference_id)
            if ref1 == ref2 and not any(h in ref1.lower() for h in ("chr", "mm39", "snapgene", "plasmid", "mouse")):
                # PASSES TIER 1 STRICT REFERENCE_UNIQUE!
                unique_counts[ref1] += 1
                is_tier1_pass = True
                sp_id, sp_name = tax_map.get(ref1, ("NO_SPECIES_RANK", "NO_SPECIES_NAME"))
                if sp_id not in species_acc_map:
                    species_acc_map[sp_id] = set()
                species_acc_map[sp_id].add(ref1)

        # STRICT ISOLATION: If passed Tier 1, STOP! Do NOT search for species-supportive!
        if is_tier1_pass:
            continue

        # TIER 2: Evaluate remaining reads for Species-Supportive assignment
        target_refs = set()
        has_host_decoy = False

        for r in alignments:
            ref_name = bamfile.get_reference_name(r.reference_id)
            if any(h in ref_name.lower() for h in ("chr", "mm39", "snapgene", "plasmid", "mouse")):
                has_host_decoy = True
                break
            # Tier 2 relaxed criteria
            tlen = abs(r.template_length)
            if r.mapping_quality >= 5 and 150 <= tlen <= 1500 and not r.is_supplementary:
                target_refs.add(ref_name)

        if has_host_decoy:
            decoy_host_count += 1
            continue

        if not target_refs:
            continue

        # Check if all target references share the SAME species_taxid
        target_sp_taxids = set()
        for tr in target_refs:
            sp_id, sp_nm = tax_map.get(tr, ("NO_SPECIES_RANK", "NO_SPECIES_NAME"))
            target_sp_taxids.add(sp_id)

        if len(target_sp_taxids) == 1 and "NO_SPECIES_RANK" not in target_sp_taxids:
            # PASSES TIER 2 SPECIES_SUPPORTIVE!
            matched_sp_id = list(target_sp_taxids)[0]
            supportive_counts[matched_sp_id] += 1
            if matched_sp_id not in species_acc_map:
                species_acc_map[matched_sp_id] = set()
            species_acc_map[matched_sp_id].update(target_refs)
        else:
            ambiguous_count += 1

    bamfile.close()

    # Generate Refined Report Rows
    report_rows = []
    all_species = set(list(species_acc_map.keys()) + [tax_map[acc][0] for acc in unique_counts.keys() if acc in tax_map])

    human_cov = args.human_coverage if args.human_coverage > 0 else 30.0
    denom = 3.1e9 / float(human_cov)

    for sp_id in sorted(all_species):
        if sp_id.startswith("NO_"):
            continue

        acc_set = species_acc_map.get(sp_id, set())
        sp_name = "NO_SPECIES_NAME"
        for acc in acc_set:
            if acc in tax_map:
                sp_name = tax_map[acc][1]
                break

        # Best supported virus accession for this species
        best_acc = "None"
        best_u_cnt = 0
        for acc in acc_set:
            u_cnt = unique_counts[acc]
            if u_cnt >= best_u_cnt:
                best_u_cnt = u_cnt
                best_acc = acc

        u_reads = unique_counts[best_acc] if best_acc != "None" else 0
        s_reads = supportive_counts[sp_id]
        total_refined = u_reads + s_reads

        if total_refined == 0:
            continue

        # Calculate Refined Copy Number (reference-discriminating copy-number lower bound)
        viral_len = ref_lengths.get(best_acc, 10000)
        norm_cov = total_refined / float(viral_len) if viral_len > 0 else 0.0
        copy_number = norm_cov * (1.0 / denom) if denom > 0 else 0.0

        status = "CONFIRMED_UNIQUE" if u_reads > 0 else "SPECIES_COLLAPSED"

        row = [
            sample_id,
            phase_str,
            proj_str,
            best_acc,
            sp_id,
            sp_name,
            str(u_reads),
            str(s_reads),
            str(total_refined),
            f"{copy_number:.6f}",
            status
        ]
        report_rows.append(row)

    header = [
        "Sample_ID", "Phase", "Project", "Virus_Accession",
        "Species_TaxID", "Species_Name", "Reference_Unique_Reads",
        "Species_Supportive_Reads", "Total_Refined_Reads",
        "Refined_Copy_Number", "Classification_Status"
    ]

    with open(out_report_tsv, "w", encoding="utf-8") as rf:
        rf.write("\t".join(header) + "\n")
        for r in report_rows:
            rf.write("\t".join(r) + "\n")

    print(f"[SUCCESS] Refinement completed for {sample_id}: {len(report_rows)} viral species records written to {out_report_tsv}")


if __name__ == "__main__":
    main()
