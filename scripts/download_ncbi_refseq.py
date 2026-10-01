#!/usr/bin/env python3
"""
download_ncbi_refseq.py
NCBI RefSeq Complete Genome Downloader for VIROnator Refinement Module

1. Reads unique viral accession IDs from Column 2 (Virus_Accession) of master_report_file.
2. Looks up species_taxid (or genus_taxid) in viral_reference_taxonomy_index.tsv.
3. Queries NCBI Entrez API for all complete nucleotide RefSeq genomes under detected species taxIDs.
4. Caches sequence files in /mnt/disks/staff/refs/ncbi_download/ (skips re-downloading existing files).
5. Sanitizes headers using scripts/rename_fasta_contigs.sh --ncbi.
6. Writes audit log to /mnt/disks/staff/SSC_hg38_refinement/logs/refinement_ncbi_fetch.log.
"""

import os
import sys
import argparse
import subprocess
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET
import time


def parse_args():
    parser = argparse.ArgumentParser(description="Download complete NCBI RefSeq viral genomes for refinement module.")
    parser.add_argument("--master-report", required=True, help="Path to master_all_cohorts_viral_report_final.tsv")
    parser.add_argument("--taxonomy-index", required=True, help="Path to viral_reference_taxonomy_index.tsv")
    parser.add_argument("--download-dir", default="/mnt/disks/staff/refs/ncbi_download", help="Directory to store downloaded FASTAs")
    parser.add_argument("--log-dir", default=os.path.join(os.environ.get("OUTPUT_DIR", "/mnt/disks/staff"), "logs"), help="Log output directory")
    parser.add_argument("--download-rank", default="species", choices=["species", "genus"], help="Taxonomic rank for Entrez query")
    parser.add_argument("--rename-script", default="scripts/rename_fasta_contigs.sh", help="Path to rename_fasta_contigs.sh script")
    return parser.parse_args()


def load_unique_report_accessions(master_report_path):
    """Extracts unique viral accession IDs from Column 2 (Virus_Accession) of master_report_file."""
    accessions = set()
    if not os.path.exists(master_report_path):
        print(f"[WARNING] Master report file not found: {master_report_path}")
        return accessions

    with open(master_report_path, "r", encoding="utf-8", errors="replace") as f:
        header = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if header is None and ("Virus_Accession" in line or "Sample_ID" in line):
                header = [p.lower() for p in parts]
                continue

            acc_col_idx = 1 # Column 2 (0-indexed 1) in master report
            if header and "virus_accession" in header:
                acc_col_idx = header.index("virus_accession")

            if len(parts) > acc_col_idx:
                acc = parts[acc_col_idx].strip()
                if acc and acc.lower() not in ("none", "0", "virus_accession", "n/a"):
                    accessions.add(acc)
                    # Also strip version dot if present
                    if "." in acc:
                        accessions.add(acc.split(".")[0])
    return accessions


def load_taxonomy_index(tax_index_path):
    """Loads taxonomy index mapping: accession -> {species_taxid, genus_taxid, species_name, genus_name}."""
    tax_map = {}
    if not os.path.exists(tax_index_path):
        print(f"[WARNING] Taxonomy index file not found: {tax_index_path}")
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
                
                # Extract indices dynamically
                sp_tax_idx = header.index("species_taxid") if header and "species_taxid" in header else 4
                gn_tax_idx = header.index("genus_taxid") if header and "genus_taxid" in header else 6
                sp_nm_idx = header.index("species_name") if header and "species_name" in header else 5
                gn_nm_idx = header.index("genus_name") if header and "genus_name" in header else 7

                info = {
                    "species_taxid": parts[sp_tax_idx] if len(parts) > sp_tax_idx else "NO_SPECIES_RANK",
                    "genus_taxid": parts[gn_tax_idx] if len(parts) > gn_tax_idx else "NO_GENUS_RANK",
                    "species_name": parts[sp_nm_idx] if len(parts) > sp_nm_idx else "NO_SPECIES_NAME",
                    "genus_name": parts[gn_nm_idx] if len(parts) > gn_nm_idx else "NO_GENUS_NAME"
                }
                tax_map[acc] = info
                tax_map[acc_ver] = info
    return tax_map


def entrez_search_and_fetch(taxid, rank_label, download_dir, log_file):
    """Queries NCBI Entrez for complete RefSeq genomes under taxid and downloads sequence FASTA."""
    out_fasta = os.path.join(download_dir, f"refinement_{rank_label}_{taxid}.fa")
    clean_fasta = os.path.join(download_dir, f"refinement_{rank_label}_{taxid}.renamed.fa")

    # CACHE CHECK: Skip download if clean file or raw file already exists and is non-empty
    if (os.path.exists(clean_fasta) and os.path.getsize(clean_fasta) > 0) or \
       (os.path.exists(out_fasta) and os.path.getsize(out_fasta) > 0):
        log_msg = f"[CACHE HIT] Sequences for {rank_label} TaxID {taxid} already present in {download_dir}. Skipping download.\n"
        print(log_msg.strip(), flush=True)
        with open(log_file, "a") as lf:
            lf.write(log_msg)
        return clean_fasta if os.path.exists(clean_fasta) else out_fasta

    print(f"[NCBI ENTREZ] Querying complete RefSeq genomes for {rank_label} TaxID {taxid} ...", flush=True)
    query = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP] AND \"complete genome\"[Title] NOT partial[Title]"
    
    esearch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={urllib.parse.quote(query)}&retmode=json&retmax=100"
    try:
        req = urllib.request.urlopen(esearch_url)
        data = req.read().decode("utf-8")
        import json
        res = json.loads(data)
        id_list = res.get("esearchresult", {}).get("idlist", [])
        
        if not id_list:
            # Fallback search without complete genome title filter
            query_fb = f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP]"
            esearch_url_fb = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={urllib.parse.quote(query_fb)}&retmode=json&retmax=50"
            req_fb = urllib.request.urlopen(esearch_url_fb)
            data_fb = req_fb.read().decode("utf-8")
            res_fb = json.loads(data_fb)
            id_list = res_fb.get("esearchresult", {}).get("idlist", [])

        if not id_list:
            log_msg = f"[NOT_FOUND] No RefSeq entries found on NCBI for {rank_label} TaxID {taxid}.\n"
            print(log_msg.strip(), flush=True)
            with open(log_file, "a") as lf:
                lf.write(log_msg)
            return None

        print(f"  Found {len(id_list)} NCBI RefSeq records for TaxID {taxid}. Downloading FASTA ...", flush=True)
        ids_str = ",".join(id_list)
        efetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids_str}&rettype=fasta&retmode=text"
        
        req_fetch = urllib.request.urlopen(efetch_url)
        fasta_data = req_fetch.read().decode("utf-8")

        with open(out_fasta, "w", encoding="utf-8") as out_f:
            out_f.write(fasta_data)

        log_msg = f"[SUCCESS] Downloaded {len(id_list)} sequences for {rank_label} TaxID {taxid} -> {out_fasta}\n"
        print(log_msg.strip(), flush=True)
        with open(log_file, "a") as lf:
            lf.write(log_msg)
            
        return out_fasta

    except Exception as e:
        log_msg = f"[ERROR] Failed to fetch TaxID {taxid} from NCBI: {e}\n"
        print(log_msg.strip(), flush=True)
        with open(log_file, "a") as lf:
            lf.write(log_msg)
        return None


def main():
    args = parse_args()

    os.makedirs(args.download_dir, exist_ok=True)
    os.makedirs(args.log_dir, exist_ok=True)

    log_file = os.path.join(args.log_dir, "refinement_ncbi_fetch.log")
    with open(log_file, "a") as lf:
        lf.write(f"\n=== NCBI RefSeq Fetch Started: {time.ctime()} ===\n")

    report_accs = load_unique_report_accessions(args.master_report)
    print(f"[INFO] Extracted {len(report_accs)} unique viral accessions from master report: {args.master_report}")

    tax_map = load_taxonomy_index(args.taxonomy_index)
    print(f"[INFO] Loaded taxonomy lookup index with {len(tax_map)} entries.")

    target_taxids = set()
    not_found_accs = []

    for acc in report_accs:
        if acc in tax_map:
            taxid_key = f"{args.download_rank}_taxid"
            t_id = tax_map[acc].get(taxid_key, "")
            if t_id and not t_id.startswith("NO_"):
                target_taxids.add(t_id)
            else:
                not_found_accs.append(acc)
        else:
            not_found_accs.append(acc)

    if not_found_accs:
        print(f"[WARNING] {len(not_found_accs)} accessions not found in taxonomy index: {', '.join(not_found_accs[:5])}...")
        with open(log_file, "a") as lf:
            lf.write(f"[NOT_FOUND] Accessions missing from taxonomy index: {', '.join(not_found_accs)}\n")

    print(f"[INFO] Querying NCBI for {len(target_taxids)} target {args.download_rank} TaxIDs ...")

    downloaded_fastas = []
    for t_id in sorted(target_taxids):
        res_fa = entrez_search_and_fetch(t_id, args.download_rank, args.download_dir, log_file)
        if res_fa:
            downloaded_fastas.append(res_fa)
        time.sleep(0.3) # Rate limit courtesy pause

    # Run header sanitization on raw downloaded fastas using rename_fasta_contigs.sh --ncbi
    if os.path.exists(args.rename_script):
        for raw_fa in downloaded_fastas:
            if not raw_fa.endswith(".renamed.fa"):
                prefix = raw_fa.rsplit(".", 1)[0]
                cmd_rename = f"bash \"{args.rename_script}\" --ncbi \"{raw_fa}\" \"{prefix}\" 2>/dev/null || true"
                subprocess.run(cmd_rename, shell=True)

    print(f"[SUCCESS] RefSeq genome fetching complete. {len(downloaded_fastas)} files ready in {args.download_dir}")


if __name__ == "__main__":
    main()
