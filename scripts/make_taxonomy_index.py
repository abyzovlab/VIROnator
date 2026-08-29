#!/usr/bin/env python3
"""
make_taxonomy_index.py
Standalone Utility Script for VIROnator

Builds a comprehensive 14-column viral taxonomy lookup index TSV from NCBI taxdump files
(nodes.dmp, names.dmp) and an accession-to-taxid mapping file. Automatically downloads
and caches NCBI taxdump files into /mnt/disks/staff/refs/ncbi_taxdump/ if not provided.

Output Schema (14 Columns):
1.  accession
2.  accession_version
3.  accession_taxid
4.  accession_name
5.  species_taxid
6.  species_name
7.  genus_taxid
8.  genus_name
9.  family_taxid
10. family_name
11. order_taxid
12. order_name
13. realm_taxid
14. realm_name
"""

import sys
import os
import argparse
import subprocess
import gzip


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build viral taxonomy lookup index from NCBI taxdump (nodes.dmp, names.dmp) and accession2taxid."
    )
    parser.add_argument(
        "-d", "--db-fasta", default="",
        help="Optional input database FASTA file (e.g. HumanViral_Reference_02-07-2022.fa) to filter taxonomy index to database accessions only."
    )
    parser.add_argument(
        "-a", "--accession-taxid", default="",
        help="Path to NCBI accession2taxid mapping file (defaults to /mnt/disks/staff/refs/ncbi_taxdump/nucl_gb.accession2taxid)"
    )
    parser.add_argument(
        "-n", "--nodes", default="",
        help="Path to NCBI nodes.dmp file (defaults to /mnt/disks/staff/refs/ncbi_taxdump/nodes.dmp)"
    )
    parser.add_argument(
        "-m", "--names", default="",
        help="Path to NCBI names.dmp file (defaults to /mnt/disks/staff/refs/ncbi_taxdump/names.dmp)"
    )
    parser.add_argument(
        "-t", "--taxdump-dir", default="/mnt/disks/staff/refs/ncbi_taxdump",
        help="Directory to store/cache downloaded NCBI taxdump files (default: /mnt/disks/staff/refs/ncbi_taxdump)"
    )
    parser.add_argument(
        "-o", "--output", default="config/db_metadata/viral_reference_taxonomy_index.tsv",
        help="Output taxonomy index TSV path (default: config/db_metadata/viral_reference_taxonomy_index.tsv)"
    )
    return parser.parse_args()


def download_ncbi_taxdump_if_needed(taxdump_dir):
    """Downloads NCBI taxdump and accession2taxid into taxdump_dir if not present."""
    os.makedirs(taxdump_dir, exist_ok=True)
    
    nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
    names_path = os.path.join(taxdump_dir, "names.dmp")
    acc_path = os.path.join(taxdump_dir, "nucl_gb.accession2taxid")

    if not (os.path.exists(nodes_path) and os.path.exists(names_path)):
        print(f"[INFO] Downloading NCBI taxdump.tar.gz to {taxdump_dir} ...", flush=True)
        tar_path = os.path.join(taxdump_dir, "taxdump.tar.gz")
        cmd_dl = f"wget -q -O \"{tar_path}\" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz && tar -xzf \"{tar_path}\" -C \"{taxdump_dir}\" nodes.dmp names.dmp"
        subprocess.run(cmd_dl, shell=True, check=True)

    if not os.path.exists(acc_path):
        acc_gz_path = os.path.join(taxdump_dir, "nucl_gb.accession2taxid.gz")
        print(f"[INFO] Downloading NCBI nucl_gb.accession2taxid.gz to {taxdump_dir} ...", flush=True)
        cmd_acc = f"wget -q -O \"{acc_gz_path}\" ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz && gunzip -f \"{acc_gz_path}\""
        subprocess.run(cmd_acc, shell=True, check=True)

    return nodes_path, names_path, acc_path


def extract_fasta_accessions(db_fasta_path):
    """Extracts accession IDs from FASTA header lines (>accession...)."""
    print(f"[INFO] Extracting target accessions from database FASTA: {db_fasta_path} ...", flush=True)
    target_accessions = set()
    open_fn = gzip.open if db_fasta_path.endswith(".gz") else open
    with open_fn(db_fasta_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                header_id = line[1:].strip().split()[0]
                acc = header_id.split(".")[0]
                target_accessions.add(acc)
                target_accessions.add(header_id)
    print(f"[INFO] Extracted {len(target_accessions)} target accession IDs from FASTA.", flush=True)
    return target_accessions


def load_nodes(nodes_filepath):
    """Loads parent taxid and rank mapping from NCBI nodes.dmp."""
    print(f"[INFO] Reading NCBI taxonomy nodes: {nodes_filepath} ...", flush=True)
    parent = {}
    rank = {}
    with open(nodes_filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) >= 3:
                taxid, par, rnk = parts[0], parts[1], parts[2]
                parent[taxid] = par
                rank[taxid] = rnk
    print(f"[INFO] Loaded {len(parent)} taxonomy nodes.", flush=True)
    return parent, rank


def load_names(names_filepath):
    """Loads scientific names from NCBI names.dmp."""
    print(f"[INFO] Reading NCBI taxonomy names: {names_filepath} ...", flush=True)
    names = {}
    with open(names_filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [x.strip() for x in line.split("|")]
            if len(parts) >= 4:
                taxid, nm, name_class = parts[0], parts[1], parts[3]
                if name_class == "scientific name":
                    names[taxid] = nm
    print(f"[INFO] Loaded {len(names)} scientific names.", flush=True)
    return names


def get_rank_taxid(taxid, target_rank, parent, rank):
    """Traverses parent taxonomy tree upwards to find target rank taxid."""
    cur = str(taxid)
    seen = set()

    while cur in parent and cur not in seen:
        seen.add(cur)
        if rank.get(cur) == target_rank:
            return cur
        if parent.get(cur) == cur:
            break
        cur = parent.get(cur, "")

    return ""


def taxid_or_missing(taxid, label):
    return taxid if taxid else f"NO_{label}_RANK"


def name_or_missing(taxid, label, names):
    return names.get(taxid, f"NO_{label}_NAME") if taxid else f"NO_{label}_NAME"


def main():
    args = parse_args()

    # Determine paths, downloading NCBI files into /mnt/disks/staff/refs/ncbi_taxdump if missing
    def_nodes, def_names, def_acc = download_ncbi_taxdump_if_needed(args.taxdump_dir)
    nodes_file = args.nodes if args.nodes else def_nodes
    names_file = args.names if args.names else def_names
    acc_file = args.accession_taxid if args.accession_taxid else def_acc

    target_accs = extract_fasta_accessions(args.db_fasta) if args.db_fasta and os.path.exists(args.db_fasta) else None

    parent, rank = load_nodes(nodes_file)
    names = load_names(names_file)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    print(f"[INFO] Processing accession-to-taxid mapping file: {acc_file} ...", flush=True)

    header = [
        "accession", "accession_version", "accession_taxid", "accession_name",
        "species_taxid", "species_name",
        "genus_taxid", "genus_name",
        "family_taxid", "family_name",
        "order_taxid", "order_name",
        "realm_taxid", "realm_name"
    ]

    count = 0
    with open(acc_file, "r", encoding="utf-8", errors="replace") as inp, \
         open(args.output, "w", encoding="utf-8") as out:
        
        out.write("\t".join(header) + "\n")

        for line in inp:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("accession"):
                continue

            parts = line.split()
            if len(parts) < 3:
                continue

            accession = parts[0]
            accession_version = parts[1]
            accession_taxid = parts[2]

            if target_accs is not None and (accession not in target_accs and accession_version not in target_accs):
                continue

            if accession_taxid == "ERROR_OR_NO_TAXID":
                out.write(
                    f"{accession}\t{accession_version}\t{accession_taxid}\tERROR_OR_NO_NAME\t"
                    "NO_SPECIES_RANK\tNO_SPECIES_NAME\t"
                    "NO_GENUS_RANK\tNO_GENUS_NAME\t"
                    "NO_FAMILY_RANK\tNO_FAMILY_NAME\t"
                    "NO_ORDER_RANK\tNO_ORDER_NAME\t"
                    "NO_REALM_RANK\tNO_REALM_NAME\n"
                )
                continue

            species = get_rank_taxid(accession_taxid, "species", parent, rank)
            genus = get_rank_taxid(accession_taxid, "genus", parent, rank)
            family = get_rank_taxid(accession_taxid, "family", parent, rank)
            order = get_rank_taxid(accession_taxid, "order", parent, rank)
            realm = get_rank_taxid(accession_taxid, "realm", parent, rank)

            out.write(
                f"{accession}\t{accession_version}\t{accession_taxid}\t{name_or_missing(accession_taxid, 'ACCESSION', names)}\t"
                f"{taxid_or_missing(species, 'SPECIES')}\t{name_or_missing(species, 'SPECIES', names)}\t"
                f"{taxid_or_missing(genus, 'GENUS')}\t{name_or_missing(genus, 'GENUS', names)}\t"
                f"{taxid_or_missing(family, 'FAMILY')}\t{name_or_missing(family, 'FAMILY', names)}\t"
                f"{taxid_or_missing(order, 'ORDER')}\t{name_or_missing(order, 'ORDER', names)}\t"
                f"{taxid_or_missing(realm, 'REALM')}\t{name_or_missing(realm, 'REALM', names)}\n"
            )
            count += 1
            if count % 100000 == 0:
                print(f"  Processed {count} entries...", flush=True)

    print(f"[SUCCESS] Taxonomy index successfully built with {count} records: {args.output}")


if __name__ == "__main__":
    main()
