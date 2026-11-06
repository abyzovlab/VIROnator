#!/usr/bin/env python3
"""
scripts/tag_kraken_fasta.py
===========================
Pre-scans composite reference FASTA files and formats sequence headers with
Kraken2 taxonomy IDs (|kraken:taxid|<TAXID>) before adding to Kraken2 library.

Taxonomy Assignment Rules:
- Human contigs (chr*, HLA-*, hs38d5*, decoy*)    -> TaxID 9606 (Homo sapiens)
- Mouse contigs (mm39*)                          -> TaxID 10090 (Mus musculus)
- Plasmid contigs (SnapGene*)                    -> TaxID 32630 (Synthetic construct)
- Viral contigs (NC_*, AF_*, etc.)               -> Mapped via NCBI accession2taxid (or fallback 10239 Viruses)

Outputs:
1. FORMATTED_FASTA with |kraken:taxid|<ID> appended to sequence IDs
2. seqid2taxid.map mapping sequence IDs to TaxIDs
"""

import os
import sys

def main():
    if len(sys.argv) < 5:
        print("Usage: tag_kraken_fasta.py <fasta_in> <fasta_out> <map_out> <tax_dir>", file=sys.stderr)
        sys.exit(1)

    fasta_in = sys.argv[1]
    fasta_out = sys.argv[2]
    map_out = sys.argv[3]
    tax_dir = sys.argv[4]

    print(f"[PRE-SCAN] Scanning target sequence headers from {fasta_in}...", flush=True)
    viral_accs = set()
    with open(fasta_in, "r") as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip().split()[0]
                if not h.startswith(("chr", "HLA-", "hs38d5", "decoy", "mm39", "SnapGene")):
                    acc = h.split(".")[0]
                    viral_accs.add(acc)
                    viral_accs.add(h)

    print(f"[PRE-SCAN] Detected {len(viral_accs)} candidate viral accessions.", flush=True)

    # Load accession to taxid from NCBI mapping files
    viral_map = {}
    for acc_file in ["nucl_gb.accession2taxid", "nucl_wgs.accession2taxid"]:
        p = os.path.join(tax_dir, acc_file)
        if os.path.exists(p):
            print(f"[PRE-SCAN] Searching {acc_file} for viral accessions...", flush=True)
            with open(p, "r") as f:
                for line in f:
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        a, v, t = parts[0], parts[1], parts[2]
                        if a in viral_accs:
                            viral_map[a] = t
                        if v in viral_accs:
                            viral_map[v] = t

    print(f"[PRE-SCAN] Successfully resolved {len(viral_map)} viral accessions to NCBI TaxIDs.", flush=True)

    # Tag each header: Human=9606, Mouse=10090, Plasmid=32630, Virus=viral_taxid (or 10239 root)
    n_human, n_mouse, n_plasmid, n_viral = 0, 0, 0, 0
    with open(fasta_in, "r") as fin, open(fasta_out, "w") as fout, open(map_out, "w") as fmap:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip()
                seqid = header.split()[0]
                rest = header[len(seqid):].strip()

                if seqid.startswith(("chr", "HLA-", "hs38d5", "decoy")):
                    taxid = "9606"   # Homo sapiens
                    n_human += 1
                elif seqid.startswith("mm39"):
                    taxid = "10090"  # Mus musculus
                    n_mouse += 1
                elif seqid.startswith("SnapGene"):
                    taxid = "32630"  # Synthetic construct / cloning vector
                    n_plasmid += 1
                else:
                    acc = seqid.split(".")[0]
                    taxid = viral_map.get(seqid, viral_map.get(acc, "10239"))  # Specific virus or root 10239
                    n_viral += 1

                new_header = f">{seqid}|kraken:taxid|{taxid} {rest}".strip()
                fout.write(new_header + "\n")
                fmap.write(f"{seqid}\t{taxid}\n")
                fmap.write(f"{seqid}|kraken:taxid|{taxid}\t{taxid}\n")
            else:
                fout.write(line)

    print(f"[PRE-SCAN] Tagging complete: {n_human} Human (9606), {n_mouse} Mouse (10090), {n_plasmid} Plasmids (32630), {n_viral} Viruses.", flush=True)


if __name__ == "__main__":
    main()
