#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  clean_reference_names_and_dedup.sh names INPUT.txt OUTPUT_PREFIX
  clean_reference_names_and_dedup.sh fasta INPUT.fa OUTPUT_PREFIX [none|exact|strand|circular]

Modes for FASTA input:
  none      Rename headers only; do not remove duplicate sequences.
  exact     Remove byte-identical sequences after whitespace/case normalization.
  strand    Remove exact duplicates and full reverse-complement duplicates.
             Recommended for a mixed double-stranded DNA reference database.
  circular  Also regard circular rotations as duplicates, in either orientation.
             Use only for a FASTA containing complete circular molecules, or after
             separating whole circular plasmids from linear fragments/constructs.

Outputs for names input:
  PREFIX.cleaned_names.txt
  PREFIX.rename_map.tsv
  PREFIX.name_collisions.tsv

Outputs for FASTA input:
  PREFIX.renamed.fa
  PREFIX.dedup.fa
  PREFIX.rename_map.tsv
  PREFIX.name_collisions.tsv
  PREFIX.duplicates.tsv
USAGE
}

if [[ $# -lt 3 || $# -gt 4 ]]; then
    usage >&2
    exit 2
fi

kind=$1
input=$2
prefix=$3
mode=${4:-strand}

if [[ ! -r "$input" ]]; then
    printf 'Error: cannot read input file: %s\n' "$input" >&2
    exit 1
fi

case "$kind" in
    names)
        [[ $# -eq 3 ]] || { usage >&2; exit 2; }
        mode=none
        ;;
    fasta)
        case "$mode" in
            none|exact|strand|circular) ;;
            *) printf 'Error: unknown deduplication mode: %s\n' "$mode" >&2; usage >&2; exit 2 ;;
        esac
        ;;
    *)
        printf 'Error: first argument must be names or fasta.\n' >&2
        usage >&2
        exit 2
        ;;
esac

python3 - "$kind" "$input" "$prefix" "$mode" <<'PY'
from __future__ import annotations

import re
import sys
import unicodedata
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Tuple, Dict, Set

kind, input_name, prefix_name, dedup_mode = sys.argv[1:5]
input_path = Path(input_name)
prefix = Path(prefix_name)
prefix.parent.mkdir(parents=True, exist_ok=True)

# Semantic replacements are performed before generic punctuation cleanup.
GREEK = str.maketrans({
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
    "ε": "epsilon", "κ": "kappa", "λ": "lambda", "μ": "mu",
    "π": "pi", "σ": "sigma", "τ": "tau", "φ": "phi", "ω": "omega",
    "Α": "Alpha", "Β": "Beta", "Γ": "Gamma", "Δ": "Delta",
    "Κ": "Kappa", "Λ": "Lambda", "Μ": "Mu", "Π": "Pi",
    "Σ": "Sigma", "Τ": "Tau", "Φ": "Phi", "Ω": "Omega",
})

DASHES = str.maketrans({
    "−": "-", "–": "-", "—": "-", "‐": "-", "‑": "-",
    "＋": "+", "＆": "&", "％": "%", "＃": "#", "＠": "@",
})


def clean_name(raw: str) -> str:
    """Return a readable identifier containing only A-Z, a-z, 0-9, '_' and '-'."""
    s = raw.strip().translate(GREEK).translate(DASHES)

    # Remove trademark notation rather than retaining it as biological metadata.
    s = re.sub(r"\((?:TM|R)\)", "", s)
    s = s.replace("™", "").replace("®", "").replace("©", "")

    # Preserve meaningful plus/minus labels before deleting parentheses.
    s = s.replace("(+)", "_plus_").replace("(-)", "_minus_")
    s = re.sub(r"\(([^()]*)\+\)", r"_\1_plus_", s)
    s = re.sub(r"\(([^()]*)-\)", r"_\1_minus_", s)

    # A comma between digits is a thousands separator, not a field boundary.
    s = re.sub(r"(?<=\d),(?=\d)", "", s)

    # Preserve prime notation in names such as 3'UTR and 5'UTR.
    s = re.sub(r"(?<=\d)'(?=[A-Za-z])", "prime", s)

    # Meaning-bearing symbols.
    s = s.replace("+", "_plus_")
    s = s.replace("&", "_and_")
    s = s.replace("%", "_percent_")
    s = s.replace("#", "_num_")
    s = s.replace("@", "_at_")
    s = s.replace("±", "_plus_minus_")

    # Transliterate accented Latin characters, then remove untransliterable marks.
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")

    # Periods, whitespace, slashes, brackets, quotes, colons, etc. become separators.
    # Existing hyphens and underscores are retained.
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", s)

    # Collapse a run of separators. A run containing an underscore becomes one
    # underscore; a run containing only hyphens becomes one hyphen.
    s = re.sub(
        r"[-_]+",
        lambda m: "_" if "_" in m.group(0) else "-",
        s,
    )
    s = s.strip("_-")

    return s or "unnamed_contig"


def assign_unique_names(raw_names: Iterable[str]):
    used = set()  # type: Set[str]
    next_suffix = defaultdict(lambda: 2)
    rows = []
    collisions = []

    for index, raw in enumerate(raw_names, start=1):
        base = clean_name(raw)
        candidate = base
        if candidate in used:
            n = next_suffix[base]
            while f"{base}_{n}" in used:
                n += 1
            candidate = f"{base}_{n}"
            next_suffix[base] = n + 1
            collisions.append((index, raw, base, candidate))
        used.add(candidate)
        rows.append((index, raw, base, candidate))
    return rows, collisions


def write_name_maps(rows, collisions) -> None:
    with open(f"{prefix}.rename_map.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("record_index\toriginal_header\tbase_clean_id\tfinal_clean_id\n")
        for index, raw, base, final in rows:
            out.write(f"{index}\t{raw}\t{base}\t{final}\n")

    with open(f"{prefix}.name_collisions.tsv", "w", encoding="utf-8", newline="\n") as out:
        out.write("record_index\toriginal_header\tbase_clean_id\tfinal_clean_id\n")
        for index, raw, base, final in collisions:
            out.write(f"{index}\t{raw}\t{base}\t{final}\n")


if kind == "names":
    with input_path.open(encoding="utf-8-sig") as handle:
        raw_names = [line.rstrip("\r\n") for line in handle if line.strip()]

    rows, collisions = assign_unique_names(raw_names)
    write_name_maps(rows, collisions)

    with open(f"{prefix}.cleaned_names.txt", "w", encoding="utf-8", newline="\n") as out:
        for _, _, _, final in rows:
            out.write(final + "\n")

    print(f"Input names: {len(rows)}")
    print(f"Sanitized-name collisions resolved with _2, _3, ...: {len(collisions)}")
    print(f"Wrote: {prefix}.cleaned_names.txt")
    print(f"Wrote: {prefix}.rename_map.tsv")
    print(f"Wrote: {prefix}.name_collisions.tsv")
    raise SystemExit(0)


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    header = None  # type: Optional[str]
    sequence_parts = []  # type: List[str]

    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"Empty FASTA header at line {line_number}")
                sequence_parts = []
            elif line.strip():
                if header is None:
                    raise ValueError(
                        f"Sequence text occurs before the first FASTA header at line {line_number}"
                    )
                sequence_parts.append("".join(line.split()))

    if header is not None:
        yield header, "".join(sequence_parts)


IUPAC = set("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN",
)


def normalize_sequence(seq: str, header: str) -> str:
    seq = seq.upper()
    invalid = sorted(set(seq) - IUPAC)
    if invalid:
        shown = " ".join(repr(ch) for ch in invalid[:10])
        raise ValueError(
            f"Invalid non-IUPAC character(s) in sequence {header!r}: {shown}"
        )
    if not seq:
        raise ValueError(f"Empty sequence for FASTA record {header!r}")
    return seq


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def minimal_rotation(seq: str) -> str:
    """Lexicographically minimal rotation using Booth's O(n) algorithm."""
    n = len(seq)
    if n < 2:
        return seq
    doubled = seq + seq
    i, j, k = 0, 1, 0
    while i < n and j < n and k < n:
        a = doubled[i + k]
        b = doubled[j + k]
        if a == b:
            k += 1
            continue
        if a > b:
            i = i + k + 1
            if i <= j:
                i = j + 1
        else:
            j = j + k + 1
            if j <= i:
                j = i + 1
        k = 0
    start = min(i, j)
    return doubled[start : start + n]


def canonical_key(seq: str, mode: str) -> str:
    if mode == "exact":
        return seq
    rc = reverse_complement(seq)
    if mode == "strand":
        return min(seq, rc)
    if mode == "circular":
        return min(minimal_rotation(seq), minimal_rotation(rc))
    raise AssertionError(mode)


def circular_offset(target: str, source: str) -> Optional[int]:
    """Return left-rotation offset if target is a rotation of source."""
    if len(target) != len(source):
        return None
    n = len(source)
    pos = (source + source).find(target, 0, 2 * n - 1)
    return None if pos < 0 else pos


def duplicate_relation(seq: str, representative: str, mode: str) -> Tuple[str, str]:
    if seq == representative:
        return "exact", "."
    if mode in {"strand", "circular"} and seq == reverse_complement(representative):
        return "reverse_complement", "."
    if mode == "circular":
        offset = circular_offset(seq, representative)
        if offset is not None:
            return "circular_rotation", str(offset)
        offset = circular_offset(seq, reverse_complement(representative))
        if offset is not None:
            return "circular_reverse_complement", str(offset)
    return "equivalent", "."


def write_record(handle, identifier: str, sequence: str, width: int = 80) -> None:
    handle.write(f">{identifier}\n")
    for start in range(0, len(sequence), width):
        handle.write(sequence[start : start + width] + "\n")


records_raw = list(read_fasta(input_path))
if not records_raw:
    raise ValueError("No FASTA records found")

rows, collisions = assign_unique_names(header for header, _ in records_raw)
write_name_maps(rows, collisions)

records = []  # type: List[Tuple[str, str, str]]
for (original_header, raw_seq), (_, _, _, clean_id) in zip(records_raw, rows):
    records.append((clean_id, original_header, normalize_sequence(raw_seq, original_header)))

with open(f"{prefix}.renamed.fa", "w", encoding="ascii", newline="\n") as out:
    for clean_id, _, seq in records:
        write_record(out, clean_id, seq)

if dedup_mode == "none":
    representatives = records
    duplicates = []  # type: List[Tuple[str, str, str, str, int]]
else:
    first_by_key = {}  # type: Dict[str, Tuple[str, str]]
    representatives = []
    duplicates = []

    for clean_id, original_header, seq in records:
        key = canonical_key(seq, dedup_mode)
        previous = first_by_key.get(key)
        if previous is None:
            first_by_key[key] = (clean_id, seq)
            representatives.append((clean_id, original_header, seq))
        else:
            representative_id, representative_seq = previous
            relation, offset = duplicate_relation(seq, representative_seq, dedup_mode)
            duplicates.append((clean_id, representative_id, relation, offset, len(seq)))

with open(f"{prefix}.dedup.fa", "w", encoding="ascii", newline="\n") as out:
    for clean_id, _, seq in representatives:
        write_record(out, clean_id, seq)

with open(f"{prefix}.duplicates.tsv", "w", encoding="utf-8", newline="\n") as out:
    out.write("duplicate_id\trepresentative_id\trelation\trotation_offset\tlength\n")
    for duplicate_id, representative_id, relation, offset, length in duplicates:
        out.write(
            f"{duplicate_id}\t{representative_id}\t{relation}\t{offset}\t{length}\n"
        )

print(f"Input FASTA records: {len(records)}")
print(f"Deduplication mode: {dedup_mode}")
print(f"Sanitized-name collisions resolved with _2, _3, ...: {len(collisions)}")
print(f"Sequence duplicates removed: {len(duplicates)}")
print(f"Records retained: {len(representatives)}")
print(f"Wrote: {prefix}.renamed.fa")
print(f"Wrote: {prefix}.dedup.fa")
print(f"Wrote: {prefix}.rename_map.tsv")
print(f"Wrote: {prefix}.name_collisions.tsv")
print(f"Wrote: {prefix}.duplicates.tsv")
PY
