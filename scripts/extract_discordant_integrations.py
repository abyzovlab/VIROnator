#!/usr/bin/env python3
"""
extract_discordant_integrations.py
VIROnator Module 10: Discordant Reads & Viral Integration Breakpoints

Performs:
1. Extraction of chimeric reads from unfiltered CRAM (clean or raw):
   - Discordant read pairs (one mate Human, one mate Viral)
   - Soft-clipped / split reads with cross-species SA alignments (Primary Human + SA Viral OR Primary Viral + SA Human)
2. Rightward sliding window clustering (Seed / Anchor approach) to group reads supporting the same physical integration locus.
3. Accurate breakpoint estimation:
   - Exact single-nucleotide junction when split reads are present.
   - Median inner/outer boundary estimation when supported by discordant pairs.
4. Nearest human gene annotation:
   - Intragenic gene name if inside, or NA with closest upstream (-) and downstream (+) gene distances.
5. Exporting:
   - exogeneSR_discordant_pairs.cram & .tsv
   - exogeneSR_softclipped_split.cram & .tsv
   - {sample}_integration_breakpoints.tsv
   - Optional Dual Genomic Arc / Locus coverage visualization plots.
"""

import sys
import os
import argparse
import subprocess
import collections
import bisect
import urllib.request
import gzip


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract discordant pairs & split reads and discover viral integration breakpoints."
    )
    # Target / Sample identifiers
    parser.add_argument("--sample-id", required=True, help="Sample identifier")
    parser.add_argument("--phase", default="", help="Phase identifier")
    parser.add_argument("--project", default="", help="Project identifier")
    
    # Input files & directories
    parser.add_argument("--input-cram", required=True, help="Path to input unfiltered CRAM (clean or raw)")
    parser.add_argument("--reference", required=True, help="Path to reference FASTA used to decode input CRAM")
    parser.add_argument("--human-contigs", required=True, help="Path to text file containing valid human contig names")
    parser.add_argument("--viral-contigs", required=True, help="Path to text file containing valid viral contig names")
    parser.add_argument("--gene-annotation", default="", help="Path to GENCODE/RefSeq genes BED file (optional)")
    parser.add_argument("--rename-map", default="", help="Path to viral renamed_map.tsv (optional)")
    
    # Output directories
    parser.add_argument("--vironator-dir", required=True, help="Sample output directory in vironator directory")
    parser.add_argument("--breakpoints-dir", default="", help="Global/cohort breakpoints output directory (optional)")
    parser.add_argument("--samtools-bin", default="samtools", help="Path to samtools binary")
    
    # Quality & Filtering thresholds
    parser.add_argument("--discordant-min-mapq", type=int, default=20, help="Min MAPQ for discordant mates")
    parser.add_argument("--discordant-min-aligned-len", type=int, default=35, help="Min aligned bases for discordant read")
    parser.add_argument("--softclip-min-primary-mapq", type=int, default=10, help="Min MAPQ for split read primary")
    parser.add_argument("--softclip-min-sa-mapq", type=int, default=10, help="Min MAPQ for split read SA alignment")
    parser.add_argument("--softclip-min-clip-len", type=int, default=25, help="Min length of soft-clipped portion")
    
    # Clustering thresholds
    parser.add_argument("--cluster-window", type=int, default=500, help="Rightward sliding window size in bp")
    parser.add_argument("--min-support-reads", type=int, default=2, help="Min supporting reads to output breakpoint")
    
    # Plotting options
    parser.add_argument("--plot-integrations", default="on", choices=["on", "off"], help="Generate integration plots")
    parser.add_argument("--plot-min-reads", type=int, default=3, help="Min reads threshold to trigger locus plot")
    parser.add_argument("--plot-top-loci", type=int, default=5, help="Max top candidate loci to plot")
    
    return parser.parse_args()


def load_contig_set(filepath):
    """Loads contig names from 1-column or tab-delimited contig file."""
    contigs = set()
    if not filepath or not os.path.exists(filepath):
        return contigs
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                contig = line.split("\t")[0].split()[0]
                contigs.add(contig)
    return contigs


def load_viral_names(rename_map_path):
    """Loads accession -> short/sanitized viral name from rename_map.tsv."""
    names = {}
    if not rename_map_path:
        return names
    # Check for _short.tsv companion first if present
    short_candidate = rename_map_path[:-4] + "_short.tsv" if rename_map_path.endswith(".tsv") else rename_map_path + "_short.tsv"
    target_path = short_candidate if os.path.exists(short_candidate) else rename_map_path
    if not os.path.exists(target_path):
        return names

    with open(target_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("final_clean_id"):
                continue
            parts = line.split("\t")
            if len(parts) >= 3 and parts[2].strip():
                names[parts[0].strip()] = parts[2].strip()
            elif len(parts) >= 2:
                clean_id = parts[0].strip()
                desc = parts[1].strip()
                if desc.startswith(clean_id):
                    desc = desc[len(clean_id):].strip()
                names[clean_id] = desc or clean_id
    return names


def ensure_gene_annotation_file(annot_path, default_dir="/mnt/disks/staff/refs"):
    """
    Checks if human gene annotation BED file exists. If missing or empty,
    attempts to download GENCODE basic genes BED into default_dir.
    """
    if annot_path and os.path.exists(annot_path):
        return annot_path
    
    candidate = os.path.join(default_dir, "gencode.v38.annotation.genes.bed")
    if os.path.exists(candidate):
        return candidate

    print(f"[INFO] Gene annotation file not found. Attempting download to {candidate} ...")
    try:
        os.makedirs(default_dir, exist_ok=True)
        url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.basic.annotation.gff3.gz"
        tmp_gz = os.path.join(default_dir, "gencode.v38.basic.annotation.gff3.gz")
        urllib.request.urlretrieve(url, tmp_gz)
        with gzip.open(tmp_gz, "rt", encoding="utf-8") as gz_in, open(candidate, "w", encoding="utf-8") as bed_out:
            for line in gz_in:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 9 and parts[2] == "gene":
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attrs = parts[8]
                    # Extract gene_name
                    gene_name = "Unknown"
                    for attr in attrs.split(";"):
                        if attr.startswith("gene_name="):
                            gene_name = attr.split("=")[1]
                            break
                    bed_out.write(f"{chrom}\t{start}\t{end}\t{gene_name}\t0\t{strand}\n")
        if os.path.exists(tmp_gz):
            os.remove(tmp_gz)
        print(f"[SUCCESS] Downloaded and processed gene annotations: {candidate}")
        return candidate
    except Exception as e:
        print(f"[WARNING] Could not auto-download gene annotation BED: {e}")
        return ""


def load_gene_intervals(annot_bed_path):
    """
    Loads gene intervals into a sorted list per chromosome for binary search:
    genes_by_chrom[chrom] = sorted list of (start, end, gene_name, strand)
    """
    genes_by_chrom = collections.defaultdict(list)
    if not annot_bed_path or not os.path.exists(annot_bed_path):
        return genes_by_chrom

    with open(annot_bed_path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    gname = parts[3]
                    strand = parts[5] if len(parts) >= 6 else "+"
                    genes_by_chrom[chrom].append((start, end, gname, strand))
                except ValueError:
                    continue

    for chrom in genes_by_chrom:
        genes_by_chrom[chrom].sort(key=lambda x: (x[0], x[1]))
    return genes_by_chrom


def annotate_genomic_point(chrom, pos, genes_by_chrom):
    """
    Finds if pos falls within a gene (intragenic), and determines nearest
    upstream (-) and downstream (+) elements with distances.
    Returns: (intragenic_element, upstream_element, upstream_distance_bp, downstream_element, downstream_distance_bp)
    """
    if chrom not in genes_by_chrom or not genes_by_chrom[chrom]:
        return ("NA", "NA", "NA", "NA", "NA")

    gene_list = genes_by_chrom[chrom]
    starts = [g[0] for g in gene_list]
    idx = bisect.bisect_right(starts, pos)

    # 1. Check for overlapping intragenic gene
    intragenic = []
    # Check candidates around idx
    check_start = max(0, idx - 10)
    check_end = min(len(gene_list), idx + 10)
    for i in range(check_start, check_end):
        g_start, g_end, g_name, g_str = gene_list[i]
        if g_start <= pos <= g_end:
            intragenic.append(g_name)

    intragenic_str = ",".join(sorted(set(intragenic))) if intragenic else "NA"

    # 2. Nearest Upstream Element (gene ending before pos)
    upstream_gene = "NA"
    upstream_dist = "NA"
    min_up_dist = float("inf")
    for i in range(idx - 1, -1, -1):
        g_start, g_end, g_name, g_str = gene_list[i]
        if g_end < pos:
            dist = pos - g_end
            if dist < min_up_dist:
                min_up_dist = dist
                upstream_gene = g_name
                upstream_dist = f"-{dist}"
            break

    # 3. Nearest Downstream Element (gene starting after pos)
    downstream_gene = "NA"
    downstream_dist = "NA"
    min_down_dist = float("inf")
    for i in range(idx, len(gene_list)):
        g_start, g_end, g_name, g_str = gene_list[i]
        if g_start > pos:
            dist = g_start - pos
            if dist < min_down_dist:
                min_down_dist = dist
                downstream_gene = g_name
                downstream_dist = f"+{dist}"
            break

    return (intragenic_str, upstream_gene, upstream_dist, downstream_gene, downstream_dist)


def parse_sa_tag(sa_str):
    """
    Parses BWA-MEM SA tag string: rname,pos,strand,CIGAR,mapQ,NM;
    Returns list of dicts.
    """
    records = []
    for entry in sa_str.strip(";").split(";"):
        if not entry:
            continue
        parts = entry.split(",")
        if len(parts) >= 5:
            records.append({
                "rname": parts[0],
                "pos": int(parts[1]),
                "strand": parts[2],
                "cigar": parts[3],
                "mapq": int(parts[4]),
            })
    return records


def get_aligned_length_from_cigar(cigar_str):
    """Calculates number of aligned bases (M, =, X) from CIGAR string."""
    import re
    aligned = 0
    for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str):
        if op in ("M", "=", "X"):
            aligned += int(length)
    return aligned


def get_soft_clipped_length(cigar_str):
    """Returns max soft-clipped length from CIGAR string."""
    import re
    clips = [int(l) for l, op in re.findall(r"(\d+)([S])", cigar_str)]
    return max(clips) if clips else 0


def determine_orientation_class(h_strand, v_strand):
    """
    Determines structural orientation class:
    - colinear: (+ , -) or (- , +)
    - inverted: (+ , +) or (- , -)
    """
    if (h_strand == "+" and v_strand == "-") or (h_strand == "-" and v_strand == "+"):
        return "colinear"
    else:
        return "inverted"


def calculate_confidence(n_discordant, n_split):
    """
    Confidence Formula:
    - High:   n_split >= 1 AND total_reads >= 3
    - Medium: (n_split >= 1 AND total_reads < 3) OR (n_split == 0 AND n_discordant >= 3)
    - Low:    n_split == 0 AND n_discordant < 3
    """
    n_total = n_discordant + n_split
    if n_split >= 1 and n_total >= 3:
        return "High"
    elif (n_split >= 1 and n_total < 3) or (n_split == 0 and n_discordant >= 3):
        return "Medium"
    else:
        return "Low"


def extract_chimeric_reads(input_cram, reference, samtools_bin, human_set, viral_set,
                           min_disc_mapq=20, min_sa_mapq=10, min_clip_len=25):
    """
    Single-pass extraction using samtools view streaming from CRAM.
    Returns:
    - discordant_pairs: list of dicts
    - split_reads: list of dicts
    - sam_header_lines: list of header lines
    """
    print(f"[INFO] Streaming reads from CRAM: {input_cram} ...")
    cmd = f"{samtools_bin} view -h --reference {reference} {input_cram}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=65536)

    sam_header_lines = []
    # Temporary buffer for mates: read_qname -> list of records
    pair_buffer = collections.defaultdict(list)
    discordant_pairs = []
    split_reads = []

    for line in proc.stdout:
        if line.startswith("@"):
            sam_header_lines.append(line)
            continue

        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 11:
            continue

        qname = parts[0]
        flag = int(parts[1])
        rname = parts[2]
        pos = int(parts[3])
        mapq = int(parts[4])
        cigar = parts[5]
        rnext = parts[6]
        pnext = int(parts[7])
        tlen = int(parts[8])
        seq = parts[9]
        qual = parts[10]
        tags = parts[11:]

        # Check for SA tag (Split reads)
        sa_str = ""
        for tag in tags:
            if tag.startswith("SA:Z:"):
                sa_str = tag[5:]
                break

        # Strand of current read
        read_strand = "-" if (flag & 0x10) else "+"

        # -------------------------------------------------------------
        # 1. Evaluate Split / Soft-clipped read with cross-species SA
        # -------------------------------------------------------------
        if sa_str and (not (flag & 0x4)): # primary/mapped
            sa_records = parse_sa_tag(sa_str)
            for sa in sa_records:
                sa_rname = sa["rname"]
                sa_mapq = sa["mapq"]
                sa_cigar = sa["cigar"]
                clip_len = get_soft_clipped_length(cigar)

                # Case 1: Primary Human, SA Viral
                is_prim_human_sa_viral = (rname in human_set and sa_rname in viral_set)
                # Case 2: Primary Viral, SA Human
                is_prim_viral_sa_human = (rname in viral_set and sa_rname in human_set)

                if (is_prim_human_sa_viral or is_prim_viral_sa_human) and mapq >= min_sa_mapq and sa_mapq >= min_sa_mapq:
                    h_chr = rname if is_prim_human_sa_viral else sa_rname
                    h_pos = pos if is_prim_human_sa_viral else sa["pos"]
                    v_chr = sa_rname if is_prim_human_sa_viral else rname
                    v_pos = sa["pos"] if is_prim_human_sa_viral else pos
                    h_strand = read_strand if is_prim_human_sa_viral else sa["strand"]
                    v_strand = sa["strand"] if is_prim_human_sa_viral else read_strand

                    split_reads.append({
                        "qname": qname,
                        "sam_line": line,
                        "human_chr": h_chr,
                        "human_pos": h_pos,
                        "human_strand": h_strand,
                        "viral_chr": v_chr,
                        "viral_pos": v_pos,
                        "viral_strand": v_strand,
                        "primary_rname": rname,
                        "primary_mapq": mapq,
                        "sa_rname": sa_rname,
                        "sa_mapq": sa_mapq,
                        "clip_len": clip_len,
                        "orientation": determine_orientation_class(h_strand, v_strand)
                    })
                    break

        # -------------------------------------------------------------
        # 2. Evaluate Discordant Read Pairs
        # -------------------------------------------------------------
        # Exclude: unmapped (0x4), secondary (0x100), supplementary (0x800), optical/PCR duplicate (0x400)
        if not (flag & (0x4 | 0x100 | 0x800 | 0x400)):
            # Check if one is human and mate is viral
            actual_rnext = rname if rnext == "=" else rnext
            is_disc = (rname in human_set and actual_rnext in viral_set) or (rname in viral_set and actual_rnext in human_set)
            if is_disc and mapq >= min_disc_mapq:
                read_record = {
                    "qname": qname,
                    "flag": flag,
                    "rname": rname,
                    "pos": pos,
                    "end": pos + max(1, get_aligned_length_from_cigar(cigar)) - 1,
                    "mapq": mapq,
                    "strand": read_strand,
                    "rnext": actual_rnext,
                    "pnext": pnext,
                    "cigar": cigar,
                    "aligned_len": get_aligned_length_from_cigar(cigar),
                    "sam_line": line
                }
                pair_buffer[qname].append(read_record)

    proc.stdout.close()
    proc.wait()

    # Form paired discordant records
    for qname, mates in pair_buffer.items():
        if len(mates) >= 2:
            m1 = mates[0]
            m2 = mates[1]
            if (m1["rname"] in human_set and m2["rname"] in viral_set):
                h_mate = m1
                v_mate = m2
            elif (m2["rname"] in human_set and m1["rname"] in viral_set):
                h_mate = m2
                v_mate = m1
            else:
                continue

            orient = determine_orientation_class(h_mate["strand"], v_mate["strand"])
            discordant_pairs.append({
                "qname": qname,
                "human_chr": h_mate["rname"],
                "human_start": h_mate["pos"],
                "human_end": h_mate["end"],
                "human_strand": h_mate["strand"],
                "human_mapq": h_mate["mapq"],
                "human_aligned_len": h_mate["aligned_len"],
                "human_sam_line": h_mate["sam_line"],
                "viral_chr": v_mate["rname"],
                "viral_start": v_mate["pos"],
                "viral_end": v_mate["end"],
                "viral_strand": v_mate["strand"],
                "viral_mapq": v_mate["mapq"],
                "viral_aligned_len": v_mate["aligned_len"],
                "viral_sam_line": v_mate["sam_line"],
                "orientation": orient
            })

    print(f"[INFO] Extracted {len(discordant_pairs)} valid discordant read pairs and {len(split_reads)} cross-species split reads.")
    return discordant_pairs, split_reads, sam_header_lines


def cluster_discordant_and_split_reads(discordant_pairs, split_reads, window_size=500):
    """
    Rightward sliding window clustering (Seed / Anchor approach):
    Groups reads that share same (human_chr, viral_chr, orientation_class)
    and fall within window_size on human and viral genomes.
    """
    # Group evidence items by (human_chr, viral_chr, orientation)
    grouped_evidence = collections.defaultdict(list)

    for p in discordant_pairs:
        key = (p["human_chr"], p["viral_chr"], p["orientation"])
        grouped_evidence[key].append({
            "type": "discordant",
            "human_start": p["human_start"],
            "human_end": p["human_end"],
            "viral_start": p["viral_start"],
            "viral_end": p["viral_end"],
            "data": p
        })

    for s in split_reads:
        key = (s["human_chr"], s["viral_chr"], s["orientation"])
        grouped_evidence[key].append({
            "type": "split",
            "human_start": s["human_pos"],
            "human_end": s["human_pos"] + 1,
            "viral_start": s["viral_pos"],
            "viral_end": s["viral_pos"] + 1,
            "data": s
        })

    clusters = []

    for (h_chr, v_chr, orient), items in grouped_evidence.items():
        # Sort left-to-right by human_start, then human_end
        items.sort(key=lambda x: (x["human_start"], x["human_end"]))

        active_clusters = []

        for item in items:
            assigned = False
            i_h_start = item["human_start"]
            i_h_end = item["human_end"]
            i_v_start = item["viral_start"]
            i_v_end = item["viral_end"]

            # Try to merge into an existing active cluster (sliding rightward)
            for c in active_clusters:
                # Human check: must not exceed human_max_end + window_size
                h_dist = i_h_start - c["human_max_end"]
                # Viral check: nearest base distance within window_size
                v_dist = max(0, max(c["viral_min_start"] - i_v_end, i_v_start - c["viral_max_end"]))

                if h_dist <= window_size and v_dist <= window_size:
                    # Append and slide window
                    c["items"].append(item)
                    c["human_max_end"] = max(c["human_max_end"], i_h_end)
                    c["viral_min_start"] = min(c["viral_min_start"], i_v_start)
                    c["viral_max_end"] = max(c["viral_max_end"], i_v_end)
                    assigned = True
                    break

            if not assigned:
                # Seed a new cluster
                new_cluster = {
                    "human_chr": h_chr,
                    "viral_chr": v_chr,
                    "orientation": orient,
                    "human_min_start": i_h_start,
                    "human_max_end": i_h_end,
                    "viral_min_start": i_v_start,
                    "viral_max_end": i_v_end,
                    "items": [item]
                }
                active_clusters.append(new_cluster)

        clusters.extend(active_clusters)

    return clusters


def write_cram(sam_lines, header_lines, out_cram, reference, samtools_bin):
    """Writes SAM lines to sorted and indexed CRAM file."""
    os.makedirs(os.path.dirname(os.path.abspath(out_cram)), exist_ok=True)
    tmp_sam = out_cram + ".tmp.sam"
    with open(tmp_sam, "w", encoding="utf-8") as f:
        for h in header_lines:
            f.write(h if h.endswith("\n") else h + "\n")
        for s in sam_lines:
            f.write(s if s.endswith("\n") else s + "\n")

    cmd_sort = f"{samtools_bin} sort -O CRAM --reference {reference} -o {out_cram} {tmp_sam}"
    cmd_index = f"{samtools_bin} index {out_cram}"
    subprocess.run(cmd_sort, shell=True, check=True)
    subprocess.run(cmd_index, shell=True, check=True)
    if os.path.exists(tmp_sam):
        os.remove(tmp_sam)


def generate_arc_plot(cluster, out_png, sample_id):
    """
    Generates a Dual Genomic Arc Plot illustrating host chromosome
    and viral genome coordinates linked by supporting reads.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
        h_chr = cluster["human_chr"]
        v_chr = cluster["viral_chr"]
        h_start = cluster["human_min_start"]
        h_end = cluster["human_max_end"]
        v_start = cluster["viral_min_start"]
        v_end = cluster["viral_max_end"]

        # Axis representation
        ax.plot([0, 10], [1, 1], color="#2C3E50", lw=4, solid_capstyle="round")
        ax.plot([0, 10], [4, 4], color="#E74C3C", lw=4, solid_capstyle="round")

        ax.text(5, 0.7, f"Host: {h_chr} ({h_start:,} - {h_end:,} bp)", ha="center", fontsize=11, weight="bold", color="#2C3E50")
        ax.text(5, 4.3, f"Viral: {v_chr} ({v_start:,} - {v_end:,} bp)", ha="center", fontsize=11, weight="bold", color="#E74C3C")

        # Draw connecting arcs for evidence
        n_items = len(cluster["items"])
        for i, item in enumerate(cluster["items"]):
            x1 = 2.0 + (i / max(1, n_items - 1)) * 6.0
            x2 = 2.0 + (i / max(1, n_items - 1)) * 6.0
            color = "#3498DB" if item["type"] == "discordant" else "#2ECC71"
            con = patches.ConnectionPatch(xyA=(x1, 1), xyB=(x2, 4), coordsA="data", coordsB="data",
                                          axesA=ax, axesB=ax, color=color, alpha=0.7, lw=2.0,
                                          arrowstyle="->", mutation_scale=12)
            ax.add_artist(con)

        ax.set_xlim(-1, 11)
        ax.set_ylim(0, 5)
        ax.axis("off")
        ax.set_title(f"Integration Breakpoint Arc Plot\nSample: {sample_id} | Support: {n_items} reads ({cluster['confidence']})", fontsize=12)

        # Legend
        p1 = patches.Patch(color="#3498DB", label=f"Discordant Pairs ({cluster['n_discordant']})")
        p2 = patches.Patch(color="#2ECC71", label=f"Split Reads ({cluster['n_split']})")
        ax.legend(handles=[p1, p2], loc="upper right", frameon=False, fontsize=9)

        plt.tight_layout()
        out_base = os.path.splitext(out_png)[0]
        out_tiff = out_base + ".tiff"
        out_png_file = out_base + ".png"
        fig.savefig(out_tiff, dpi=300, bbox_inches="tight")
        fig.savefig(out_png_file, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[SUCCESS] Saved Dual Genomic Arc Plot: {out_tiff} (300 DPI TIFF) and {out_png_file}")
    except Exception as e:
        print(f"[WARNING] Could not generate arc plot: {e}")


def main():
    args = parse_args()

    sample_id = args.sample_id
    phase = args.phase
    project = args.project
    human_contigs_file = args.human_contigs
    viral_contigs_file = args.viral_contigs
    input_cram = args.input_cram
    reference = args.reference
    samtools_bin = args.samtools_bin
    vironator_dir = args.vironator_dir
    breakpoints_dir = args.breakpoints_dir

    os.makedirs(vironator_dir, exist_ok=True)
    if breakpoints_dir:
        os.makedirs(breakpoints_dir, exist_ok=True)

    print("==================================================================")
    print("VIROnator Module 10: Discordant Reads & Viral Integration Breakpoints")
    print(f"Sample:        {sample_id}")
    print(f"Phase:         {phase}")
    print(f"Project:       {project}")
    print(f"Input CRAM:    {input_cram}")
    print(f"Vironator Dir: {vironator_dir}")
    print(f"Breakpoints:   {breakpoints_dir}")
    print("==================================================================")

    # 1. Load contig sets
    human_set = load_contig_set(human_contigs_file)
    viral_set = load_contig_set(viral_contigs_file)
    viral_names = load_viral_names(args.rename_map)

    print(f"[INFO] Loaded {len(human_set)} human contigs and {len(viral_set)} viral contigs.")

    # 2. Extract chimeric discordant pairs and split reads
    disc_pairs, split_reads, sam_headers = extract_chimeric_reads(
        input_cram=input_cram,
        reference=reference,
        samtools_bin=samtools_bin,
        human_set=human_set,
        viral_set=viral_set,
        min_disc_mapq=args.discordant_min_mapq,
        min_sa_mapq=args.softclip_min_sa_mapq,
        min_clip_len=args.softclip_min_clip_len
    )

    # 3. Write individual evidence CRAM & TSV files into vironator_dir
    # A. Discordant Pairs
    disc_cram = os.path.join(vironator_dir, "exogeneSR_discordant_pairs.cram")
    disc_tsv = os.path.join(vironator_dir, "exogeneSR_discordant_pairs.tsv")
    disc_sam_lines = []
    for p in disc_pairs:
        disc_sam_lines.append(p["human_sam_line"])
        disc_sam_lines.append(p["viral_sam_line"])

    if disc_sam_lines:
        write_cram(disc_sam_lines, sam_headers, disc_cram, reference, samtools_bin)
        print(f"[SUCCESS] Saved discordant pairs CRAM: {disc_cram}")

    with open(disc_tsv, "w", encoding="utf-8") as out_d:
        out_d.write("qname\thuman_chr\thuman_start\thuman_end\thuman_strand\thuman_mapq\tviral_chr\tviral_start\tviral_end\tviral_strand\tviral_mapq\torientation\n")
        for p in disc_pairs:
            out_d.write(f"{p['qname']}\t{p['human_chr']}\t{p['human_start']}\t{p['human_end']}\t{p['human_strand']}\t{p['human_mapq']}\t{p['viral_chr']}\t{p['viral_start']}\t{p['viral_end']}\t{p['viral_strand']}\t{p['viral_mapq']}\t{p['orientation']}\n")
    print(f"[SUCCESS] Saved discordant pairs TSV: {disc_tsv}")

    # B. Split / Soft-clipped reads
    split_cram = os.path.join(vironator_dir, "exogeneSR_softclipped_split.cram")
    split_tsv = os.path.join(vironator_dir, "exogeneSR_softclipped_split.tsv")
    split_sam_lines = [s["sam_line"] for s in split_reads]

    if split_sam_lines:
        write_cram(split_sam_lines, sam_headers, split_cram, reference, samtools_bin)
        print(f"[SUCCESS] Saved soft-clipped split reads CRAM: {split_cram}")

    with open(split_tsv, "w", encoding="utf-8") as out_s:
        out_s.write("qname\tprimary_rname\tprimary_mapq\tsa_rname\tsa_mapq\thuman_chr\thuman_pos\thuman_strand\tviral_chr\tviral_pos\tviral_strand\tclip_len\torientation\n")
        for s in split_reads:
            out_s.write(f"{s['qname']}\t{s['primary_rname']}\t{s['primary_mapq']}\t{s['sa_rname']}\t{s['sa_mapq']}\t{s['human_chr']}\t{s['human_pos']}\t{s['human_strand']}\t{s['viral_chr']}\t{s['viral_pos']}\t{s['viral_strand']}\t{s['clip_len']}\t{s['orientation']}\n")
    print(f"[SUCCESS] Saved soft-clipped split reads TSV: {split_tsv}")

    # 4. Perform Rightward Sliding Window Clustering
    clusters = cluster_discordant_and_split_reads(disc_pairs, split_reads, window_size=args.cluster_window)
    print(f"[INFO] Discovered {len(clusters)} raw integration clusters.")

    # 5. Check & load gene annotations
    annot_file = ensure_gene_annotation_file(args.gene_annotation)
    gene_intervals = load_gene_intervals(annot_file)
    if gene_intervals:
        print(f"[INFO] Loaded gene annotations for {len(gene_intervals)} chromosomes from {annot_file}.")

    # 6. Process clusters and calculate breakpoints, confidence, annotations
    breakpoint_records = []

    for c in clusters:
        n_disc = sum(1 for it in c["items"] if it["type"] == "discordant")
        n_spl = sum(1 for it in c["items"] if it["type"] == "split")
        n_tot = n_disc + n_spl

        c["n_discordant"] = n_disc
        c["n_split"] = n_spl
        c["n_total"] = n_tot
        c["confidence"] = calculate_confidence(n_disc, n_spl)

        # Filter by min_support_reads
        if n_tot < args.min_support_reads:
            continue

        # Breakpoint coordinate calculation
        # If split reads exist: mode/exact junction
        if n_spl > 0:
            spl_human_coords = [it["data"]["human_pos"] for it in c["items"] if it["type"] == "split"]
            spl_viral_coords = [it["data"]["viral_pos"] for it in c["items"] if it["type"] == "split"]
            h_bp_est = collections.Counter(spl_human_coords).most_common(1)[0][0]
            v_bp_est = collections.Counter(spl_viral_coords).most_common(1)[0][0]
        else:
            # Midpoint / median boundary of cluster interval
            h_bp_est = (c["human_min_start"] + c["human_max_end"]) // 2
            v_bp_est = (c["viral_min_start"] + c["viral_max_end"]) // 2

        c["human_bp_est"] = h_bp_est
        c["viral_bp_est"] = v_bp_est

        # Nearest gene annotation
        intra, up_gene, up_dist, down_gene, down_dist = annotate_genomic_point(c["human_chr"], h_bp_est, gene_intervals)

        v_name = viral_names.get(c["viral_chr"], c["viral_chr"])
        h_win = f"{c['human_chr']}:{c['human_min_start']}-{c['human_max_end']}"
        v_win = f"{c['viral_chr']}:{c['viral_min_start']}-{c['viral_max_end']}"

        breakpoint_records.append({
            "sample_id": sample_id,
            "phase": phase or "unknown",
            "project": project or "base",
            "human_chr": c["human_chr"],
            "human_coord_window": h_win,
            "human_breakpoint_est": h_bp_est,
            "viral_accession": c["viral_chr"],
            "viral_name_sanitized": v_name,
            "viral_coord_window": v_win,
            "viral_breakpoint_est": v_bp_est,
            "orientation": c["orientation"],
            "discordant_read_count": n_disc,
            "split_read_count": n_spl,
            "total_evidence_reads": n_tot,
            "confidence": c["confidence"],
            "intragenic_element": intra,
            "upstream_element": up_gene,
            "upstream_distance_bp": up_dist,
            "downstream_element": down_gene,
            "downstream_distance_bp": down_dist,
            "cluster_obj": c
        })

    # Sort breakpoint records by total_evidence_reads descending
    breakpoint_records.sort(key=lambda x: x["total_evidence_reads"], reverse=True)

    # 7. Write Breakpoint TSV Report to sample vironator_dir (and optional breakpoints_dir)
    sample_bp_tsv = os.path.join(vironator_dir, f"{sample_id}_integration_breakpoints.tsv")
    target_tsvs = [sample_bp_tsv]
    if breakpoints_dir:
        target_tsvs.append(os.path.join(breakpoints_dir, f"{sample_id}_integration_breakpoints.tsv"))

    header = [
        "sample_id", "phase", "project", "human_chr", "human_coord_window", "human_breakpoint_est",
        "viral_accession", "viral_name_sanitized", "viral_coord_window", "viral_breakpoint_est",
        "orientation", "discordant_read_count", "split_read_count", "total_evidence_reads",
        "confidence", "intragenic_element", "upstream_element", "upstream_distance_bp",
        "downstream_element", "downstream_distance_bp"
    ]

    for out_tsv in target_tsvs:
        with open(out_tsv, "w", encoding="utf-8") as f_bp:
            f_bp.write("\t".join(header) + "\n")
            for r in breakpoint_records:
                row_vals = [str(r[k]) for k in header]
                f_bp.write("\t".join(row_vals) + "\n")
        print(f"[SUCCESS] Written Breakpoints TSV report: {out_tsv} ({len(breakpoint_records)} candidate loci)")

    # 8. Generate Locus Arc Plots for Top Candidate Breakpoints
    if args.plot_integrations == "on":
        top_candidates = [r for r in breakpoint_records if r["total_evidence_reads"] >= args.plot_min_reads][:args.plot_top_loci]
        print(f"[INFO] Generating integration plots for top {len(top_candidates)} candidate loci...")
        for rank, cand in enumerate(top_candidates, 1):
            plot_base = f"{sample_id}_locus_rank{rank}_{cand['human_chr']}_{cand['viral_accession']}_arc"
            plot_path_sample = os.path.join(vironator_dir, plot_base + ".png")
            generate_arc_plot(cand["cluster_obj"], plot_path_sample, sample_id)
            if breakpoints_dir:
                import shutil
                for ext in [".tiff", ".png"]:
                    src = os.path.join(vironator_dir, plot_base + ext)
                    dst = os.path.join(breakpoints_dir, plot_base + ext)
                    if os.path.exists(src):
                        shutil.copy2(src, dst)

    print("\n[COMPLETED] Module 10 integration discovery finished successfully.")


if __name__ == "__main__":
    main()
