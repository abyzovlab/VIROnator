#!/usr/bin/env bash
set -euo pipefail

contig_prefix=""
is_ncbi="false"
input=""
prefix=""

while [[ "$#" -gt 0 ]]; do case $1 in
    -p|--prefix) contig_prefix="$2"; shift;;
    -n|--ncbi) is_ncbi="true";;
    *) if [[ -z "$input" ]]; then input="$1"; elif [[ -z "$prefix" ]]; then prefix="$1"; else echo "Unknown parameter: $1" >&2; exit 1; fi;;
esac; shift; done

if [[ -z "$input" || -z "$prefix" ]]; then
    echo "Usage: $0 [-p CONTIG_PREFIX] [--ncbi] INPUT.fa OUTPUT_PREFIX" >&2
    echo "Example for NCBI reference: $0 --ncbi ncbi_viral_raw.fa ncbi_viral_modified" >&2
    echo "Example for mouse reference: $0 -p \"mm39_\" mm39_raw.fa mm39_ms_modified" >&2
    exit 2
fi

[[ -r "$input" ]] || { echo "Error: cannot read $input" >&2; exit 1; }

python3 - "$input" "$prefix" "$contig_prefix" "$is_ncbi" <<'PY'
import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path

input_path = Path(sys.argv[1])
prefix = Path(sys.argv[2])
contig_prefix = sys.argv[3] if len(sys.argv) > 3 else ""
is_ncbi = (sys.argv[4].lower() == "true") if len(sys.argv) > 4 else False
prefix.parent.mkdir(parents=True, exist_ok=True)

GREEK = str.maketrans({
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
    "κ": "kappa", "λ": "lambda", "μ": "mu", "π": "pi",
    "σ": "sigma", "τ": "tau", "φ": "phi", "ω": "omega",
})
DASHES = str.maketrans({"−": "-", "–": "-", "—": "-", "‐": "-", "‑": "-"})


def clean_name(raw: str) -> str:
    s = raw.strip().translate(GREEK).translate(DASHES)
    s = re.sub(r"\((?:TM|R)\)", "", s, flags=re.IGNORECASE)
    s = s.replace("™", "").replace("®", "").replace("©", "")
    s = s.replace("(+)", "_plus_").replace("(-)", "_minus_")
    s = re.sub(r"\(([^()]*)\+\)", r"_\1_plus_", s)
    s = re.sub(r"\(([^()]*)-\)", r"_\1_minus_", s)
    s = re.sub(r"(?<=\d),(?=\d)", "", s)
    s = re.sub(r"(?<=\d)'(?=[A-Za-z])", "prime", s)
    s = s.replace("+", "_plus_").replace("&", "_and_")
    s = s.replace("%", "_percent_").replace("#", "_num_")
    s = s.replace("@", "_at_").replace("±", "_plus_minus_")
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", s)
    s = re.sub(r"[-_]+", lambda m: "_" if "_" in m.group(0) else "-", s)
    s = s.strip("_-")
    clean = s or "unnamed_contig"
    if contig_prefix and not clean.startswith(contig_prefix):
        clean = f"{contig_prefix}{clean}"
    return clean

used = set()
next_suffix = defaultdict(lambda: 2)
rows = []
collisions = []
record_index = 0

out_fasta = Path(f"{prefix}.renamed.fa")
out_map = Path(f"{prefix}.rename_map.tsv")
out_collisions = Path(f"{prefix}.name_collisions.tsv")
out_contigs = Path(f"{prefix}.contigs.txt")
out_bed = Path(f"{prefix}.bed")

contig_lengths = defaultdict(int)
current_contig = None

with input_path.open("r", encoding="utf-8-sig", newline="") as src, \
     out_fasta.open("w", encoding="utf-8", newline="") as dst:
    for line_number, line in enumerate(src, start=1):
        if not line.startswith(">"):
            dst.write(line)  # sequence lines are copied unchanged
            if current_contig:
                contig_lengths[current_contig] += len(line.strip())
            continue

        record_index += 1
        ending = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
        raw = line[1:].rstrip("\r\n")
        if not raw.strip():
            raise SystemExit(f"Error: empty FASTA header at line {line_number}")

        ncbi_desc = ""
        if is_ncbi:
            parts = raw.strip().split(maxsplit=1)
            target_raw = parts[0]
            ncbi_desc = parts[1] if len(parts) > 1 else ""
        else:
            target_raw = raw

        base = clean_name(target_raw)
        final = base
        if final in used:
            n = next_suffix[base]
            while f"{base}_{n}" in used:
                n += 1
            final = f"{base}_{n}"
            next_suffix[base] = n + 1
            collisions.append((record_index, raw, base, final, ncbi_desc))

        used.add(final)
        current_contig = final
        rows.append((final, raw))
        if final != base:
            collisions.append((final, raw))
        dst.write(">" + final + ending)

if record_index == 0:
    raise SystemExit("Error: no FASTA headers beginning with '>' were found")

with out_map.open("w", encoding="utf-8", newline="\n") as out:
    out.write("final_clean_id\toriginal_header\n")
    for row in rows:
        out.write(f"{row[0]}\t{row[1]}\n")

with out_collisions.open("w", encoding="utf-8", newline="\n") as out:
    out.write("final_clean_id\toriginal_header\n")
    for row in collisions:
        out.write(f"{row[0]}\t{row[1]}\n")

with out_contigs.open("w", encoding="utf-8", newline="\n") as out:
    for row in rows:
        out.write(str(row[0]) + "\n")

with out_bed.open("w", encoding="utf-8", newline="\n") as out:
    for row in rows:
        c_id = row[0]
        c_len = contig_lengths.get(c_id, 0)
        out.write(f"{c_id}\t0\t{c_len}\n")

print(f"Renamed FASTA: {out_fasta}")
print(f"Name map: {out_map}")
print(f"Collisions: {out_collisions} ({len(collisions)})")
print(f"Contigs List: {out_contigs}")
print(f"BED File: {out_bed}")
print(f"Records: {record_index}")
PY
